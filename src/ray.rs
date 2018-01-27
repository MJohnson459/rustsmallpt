#![cfg_attr(feature = "unstable", feature(test))]

use std::f64::consts::PI;
use rand::ThreadRng;
use rand::Rng;
use na::{Vector3};

use scene::Scene;
use sphere::Sphere;
use material::ReflectType;

// ----------- Ray --------------

#[derive(Copy, Clone, Debug)]
pub struct Ray {
    pub origin: Vector3<f64>,
    pub direction: Vector3<f64>,
}

pub fn radiance(scene: &Scene, ray: &Ray, mut depth: i32, rng: &mut ThreadRng, emission: f64) -> Vector3<f64> {
    let (intersects, closest_intersect_distance, id) = scene.intersect(ray);
    if !intersects {  // if miss, return black
        return Vector3::<f64>::zeros();
    }

    let obj: &Sphere = &scene.spheres[id];  // the hit object

    let intersect_point = ray.origin + ray.direction * closest_intersect_distance; // point on sphere where intersects
    let surface_normal = (&intersect_point - &obj.position).normalize(); // surface light_normal of intersection point
    let light_normal: Vector3<f64>; // corrected light_normal (ie internal or external intersection)
    if surface_normal.dot(&ray.direction) < 0.0 { // dot product negative if ray is internal
        light_normal = surface_normal;
    } else {
        light_normal = surface_normal * -1.0;
    }

    let mut f: Vector3<f64> = obj.color.clone();

    // Russian Roulette
    // Use maximum component (r,g,b) of the surface color
    let p: f64; // brightest colour
    if f.x > f.y && f.x > f.z {
        p = f.x;
    } else if f.y > f.z {
        p = f.y;
    } else {
        p = f.z;
    }

    // Don't do Russian Roulette until after depth 5
    depth = depth + 1;
    if depth > 5 || p == 0.0 {
        if rng.next_f64() < p {
            f = f * (1.0 / p); // normalize colour [0,1]
        } else {
            // This might be at the wrong if statement
            return obj.emission * (emission as f64);
        }
    }

    match obj.reflection {
        ReflectType::DIFF => {
            let random_angle: f64 = 2.0 * PI * rng.next_f64();
            let random: f64 = rng.next_f64();
            let r2s = random.sqrt();

            let w: Vector3<f64> = light_normal.clone();
            let u: Vector3<f64>;
            if w.x.abs() > 0.1 {
                u = (Vector3::<f64>::new(0.0, 1.0, 0.0).cross(&w)).normalize();
            } else {
                u = (Vector3::<f64>::new(1.0, 0.0, 0.0).cross(&w)).normalize();
            }
            assert!((u.length() - 1.0).abs() < 0.001);

            let v: Vector3<f64> = w.cross(&u);

            let random_direction = (u * random_angle.cos() * r2s + v * random_angle.sin() * r2s + w * (1.0 - random).sqrt()).normalize();
            assert!((random_direction.length() - 1.0).abs() < 0.001);


            let e = Vector3::<f64>::zeros();

            let fmu = f * (radiance(&scene, &Ray{origin: intersect_point, direction: random_direction}, depth, rng, 1.0));
            return obj.emission*(emission as f64) + e + fmu;
        },
        ReflectType::SPEC => {
            //println!("obj.emission SPEC, depth: {}", depth);
            let d: Vector3<f64> = ray.direction-(surface_normal*2.0*surface_normal.dot(&ray.direction));
            return obj.emission + f * (radiance(&scene, &Ray{origin: intersect_point, direction: d}, depth, rng, 1.0));
        },
            _ => {}

    }

    let refl_ray: Ray = Ray{origin: intersect_point, direction: ray.direction-surface_normal*2.0*surface_normal.dot(&ray.direction)}; // Ideal dielectric REFRACTION
    let into: bool = surface_normal.dot(&light_normal) > 0.0;

    let nc: f64 = 1.0;
    let nt: f64 = 1.5;
    let nnt: f64;
    if into { // Ray from outside going in?
        nnt = nc/nt;
    } else {
        nnt = nt/nc;
    }
    let ddn: f64 =ray.direction.dot(&light_normal);
    let cos2t: f64 = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
    if cos2t < 0.0 { // Total internal reflection
        //println!("obj.emission cos2t < 0.0");
        return obj.emission + f * (radiance(&scene, &refl_ray, depth, rng, 1.0));
    }

    let tdir: Vector3<f64>;

    if into {
        tdir = (ray.direction*nnt - surface_normal*(1.0*(ddn*nnt+cos2t.sqrt()))).normalize();
    } else {
        tdir = (ray.direction*nnt - surface_normal*(-1.0*(ddn*nnt+cos2t.sqrt()))).normalize();
    }

    let a: f64 = nt- nc;
    let b: f64 = nt + nc;
    let r0: f64 = a*a/(b*b);
    let c: f64;

    if into {
        c = 1.0 + ddn;
    } else {
        c = 1.0 - tdir.dot(&surface_normal);
    }

    let re: f64 = r0+(1.0-r0)*c.powi(5);
    let tr: f64 = 1.0 - re;
    let p: f64 = 0.25+0.5*re;
    let rp: f64 = re/p;
    let tp: f64 = tr/(1.0-p);


    if depth > 2 {
        if rng.next_f64() < p {
            obj.emission + f * (radiance(&scene, &refl_ray, depth, rng, 1.0)*rp)
        } else {
            obj.emission + f * (radiance(&scene, &Ray{origin: intersect_point, direction: tdir}, depth, rng, 1.0)*tp)
        }
    } else {
        obj.emission + f * (radiance(&scene, &refl_ray, depth, rng, 1.0)*re+radiance(&scene, &Ray{origin: intersect_point, direction: tdir}, depth, rng, 1.0)*tr)
    }
}
