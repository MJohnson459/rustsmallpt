use std::f64::consts::PI;
use rand::ThreadRng;
use rand::Rng;

use vector_3d::Vec3d;
use scene::Scene;
use sphere::Sphere;
use material::ReflectType;

// ----------- Ray --------------

#[derive(Copy, Clone, Debug)]
pub struct Ray {
    pub origin: Vec3d,
    pub direction: Vec3d,
}

fn get_brightest_color(color: &Vec3d) -> f64 {
    if color.x > color.y && color.x > color.z {
        color.x
    } else if color.y > color.z {
        color.y
    } else {
        color.z
    }
}

fn get_random_direction(rng: &mut ThreadRng, light_normal: &Vec3d) -> Vec3d {
    let random_angle = 2.0 * PI * rng.next_f64();
    let random = rng.next_f64();
    let r2s = random.sqrt();

    let w = *light_normal;
    let u =
        if w.x.abs() > 0.1 {
            (Vec3d{x:0.0, y:1.0, z:0.0}.cross(w)).normalise()
        } else {
            (Vec3d{x:1.0, y:0.0, z:0.0}.cross(w)).normalise()
        };

    let v = w.cross(u);

    (u * random_angle.cos() * r2s + v * random_angle.sin() * r2s + w * (1.0 - random).sqrt()).normalise()
}

fn get_reflected_ray(ray: &Ray, surface_normal: &Vec3d, surface_point: &Vec3d) -> Ray {
    let dir = ray.direction - (*surface_normal * 2.0 * surface_normal.dot(&ray.direction));
    Ray {
        origin: *surface_point,
        direction: dir,
    }
}

pub fn radiance(scene: &Scene, ray: &Ray, depth: i32, rng: &mut ThreadRng, emission: f64) -> Vec3d {
    let (intersects, closest_intersect_distance, id) = scene.intersect(ray);
    if !intersects {  // if miss, return black
        return Vec3d{x:0.0, y:0.0, z:0.0};
    }

    let obj: &Sphere = &scene.spheres[id];  // the hit object

    let intersect_point = ray.origin + ray.direction * closest_intersect_distance; // point on sphere where intersects
    let surface_normal = (&intersect_point - &obj.position).normalise(); // surface light_normal of intersection point

    // corrected light_normal (ie internal or external intersection)
    let light_normal =
        if surface_normal.dot(&ray.direction) < 0.0 { // dot product negative if ray is internal
            surface_normal
        } else {
            surface_normal * -1.0
        };


    // Russian Roulette
    // Use maximum component (r,g,b) of the surface color
    let brightest_color = get_brightest_color(&obj.color); // brightest colour

    // Don't do Russian Roulette until after depth 5
    if (depth > 5 || brightest_color == 0.0) && rng.next_f64() >= brightest_color {
        return obj.emission * emission;
    }

    let color = obj.color;

    match obj.reflection {
        ReflectType::DIFF => {
            let random_direction = get_random_direction(rng, &light_normal);
            let e = Vec3d::zeros();

            let fmu = color * radiance(&scene, &Ray{origin: intersect_point, direction: random_direction}, depth + 1, rng, 1.0);
            return obj.emission * emission + e + fmu;
        },
        ReflectType::SPEC => {
            return obj.emission + color * radiance(&scene, &get_reflected_ray(&ray, &surface_normal, &intersect_point), depth + 1, rng, 1.0);
        },
        ReflectType::REFR => {
            let refl_ray = get_reflected_ray(&ray, &surface_normal, &intersect_point);
            let into = surface_normal.dot(&light_normal) > 0.0;

            let air: f64 = 1.0;
            let glass: f64 = 1.5;
            let refraction = if into { air / glass } else { glass / air };

            let angle = ray.direction.dot(&light_normal);
            let cos2t = 1.0 - refraction * refraction * (1.0 - angle * angle);
            if cos2t < 0.0 {
                // Total internal reflection so all light is reflected (internally)
                return obj.emission + color * radiance(&scene, &refl_ray, depth + 1, rng, 1.0);
            }

            let refract_dir = (ray.direction * refraction - surface_normal * (if into {1.0} else {-1.0} * (angle * refraction + cos2t.sqrt()))).normalise();

            // Fresnel reflectance
            let a = glass - air;
            let b = glass + air;
            let normal_reflected = a * a / (b * b);
            let c =
                if into {
                    1.0 + angle
                } else {
                    1.0 - refract_dir.dot(&surface_normal)
                };

            let total_reflected = normal_reflected + (1.0 - normal_reflected) * c.powi(5);
            let total_refracted = 1.0 - total_reflected;
            let reflect_probability = 0.25 + 0.5 * total_reflected;

            // Weight results based on probability
            let reflect_weight = total_reflected / reflect_probability;
            let refract_weight = total_refracted / (1.0 - reflect_probability);

            // if depth is shallow, make 2 recursive calls
            if depth > 2 {
                if rng.next_f64() < reflect_probability {
                    obj.emission + color * radiance(&scene, &refl_ray, depth + 1, rng, 1.0) * reflect_weight
                } else {
                    obj.emission + color * radiance(&scene, &Ray{origin: intersect_point, direction: refract_dir}, depth + 1, rng, 1.0) * refract_weight
                }
            } else {
                obj.emission + color * radiance(&scene, &refl_ray, depth + 1, rng, 1.0) * total_reflected
                    + radiance(&scene, &Ray{origin: intersect_point, direction: refract_dir}, depth + 1, rng, 1.0) * total_refracted
            }
        },
    }
}

#[cfg(all(feature = "unstable", test))]
mod bench {
    extern crate test;
    extern crate rand;

    use super::*;
    use self::test::Bencher;

    #[bench]
    fn bench_radiance(b: &mut Bencher) {
        let scene = Scene::new2();
        let mut rng = rand::thread_rng();
        let ray = Ray{origin: Vec3d{x: 50.0, y: 50.0, z: 100.0}, direction: Vec3d{x: 0.0, y: -0.042612, z: -1.0}.normalise()};

        b.iter(|| {
            radiance(&scene, &ray, 0, &mut rng, 1.0)
        });
    }

    #[bench]
    fn bench_radiance_black(b: &mut Bencher) {
        let scene = Scene::new2();
        let mut rng = rand::thread_rng();
        let ray = Ray{origin: Vec3d{x: 50.0, y: 50.0, z: 260.0}, direction: Vec3d{x: 0.0, y: -0.042612, z: -1.0}.normalise()};

        b.iter(|| {
            radiance(&scene, &ray, 0, &mut rng, 1.0)
        });
    }
}
