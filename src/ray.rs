use std::f64::consts::PI;
use rand::ThreadRng;
use rand::Rng;
use num_traits::identities::Zero;
use cgmath::{InnerSpace, ElementWise};

use config::Config;
use material::ReflectType;
use sphere::Sphere;
use vector_3d::Vec3d;

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

fn get_random_direction(rng: &mut ThreadRng, oriented_surface_normal: &Vec3d) -> Vec3d {
    let random_angle = 2.0 * PI * rng.next_f64();
    let random = rng.next_f64();
    let r2s = random.sqrt();

    let w = *oriented_surface_normal;
    let u = if w.x.abs() > 0.1 {
        (Vec3d::new( 0.0, 1.0, 0.0).cross(w))
            .normalize()
    } else {
        (Vec3d::new( 1.0, 0.0, 0.0).cross(w))
            .normalize()
    };

    let v = w.cross(u);

    (u * random_angle.cos() * r2s + v * random_angle.sin() * r2s + w * (1.0 - random).sqrt())
        .normalize()
}

fn get_reflected_ray(ray: &Ray, surface_normal: &Vec3d, surface_point: &Vec3d) -> Ray {
    let dir = ray.direction - (*surface_normal * 2.0 * surface_normal.dot(ray.direction));
    Ray {
        origin: *surface_point,
        direction: dir,
    }
}

fn diff_radiance(
    config: &Config,
    _: &Ray,
    depth: i32,
    rng: &mut ThreadRng,
    emit: bool,
    surface_normal: &Vec3d,
    intersect_point: &Vec3d,
    obj: &Sphere,
) -> Vec3d {
    let mut e = Vec3d::new(0.0, 0.0, 0.0);
    if config.explicit_light_sampling {
        // Loop over any lights
        config
            .scene
            .spheres
            .iter()
            .enumerate()
            .for_each(|(i, sphere)| {
                if sphere.emission.x <= 0.0 && sphere.emission.y <= 0.0 && sphere.emission.z <= 0.0
                {
                    // skip non-lights
                } else {
                    // Create random direction towards sphere
                    let sw = sphere.position - *intersect_point;

                    if surface_normal.dot(sw) >= 0.0 {
                        let su = ((if sw.x.abs() > 0.1 {
                            Vec3d::new( 0.0, 1.0, 0.0)
                        } else {
                            Vec3d::new( 1.0, 0.0, 0.0)
                        }).cross(sw))
                            .normalize();
                        let sv = sw.cross(su);

                        let relative_intersect = *intersect_point - sphere.position;
                        let cos_a_max = (1.0
                            - sphere.radius * sphere.radius
                                / relative_intersect.dot(relative_intersect))
                            .sqrt();

                        let eps1 = rng.next_f64();
                        let eps2 = rng.next_f64();
                        let cos_a = 1.0 - eps1 + eps1 * cos_a_max;
                        let sin_a = (1.0 - cos_a * cos_a).sqrt();
                        let phi = 2.0 * PI * eps2;
                        let l = (su * phi.cos() * sin_a + sv * phi.sin() * sin_a + sw * cos_a)
                            .normalize();

                        // Shoot shadow ray
                        config.scene.intersect(&Ray {origin: *intersect_point, direction: l})
                            .and_then(| (_, id) | if id == i {
                                // shadow ray
                                let omega = 2.0 * PI * (1.0 - cos_a_max); // 1 / probability with respect to solid angle
                                e = e
                                    + obj.color.mul_element_wise(sphere.emission * l.dot(*surface_normal) * omega)
                                        * (1.0 / PI); // 1/pi for brdf
                                                      //println!("e: {:?}", e);
                                Some(e)
                            } else {
                                None
                            });
                    }
                }
            });
    }

    let new_emit = !config.explicit_light_sampling || e.is_zero();

    let random_direction = get_random_direction(rng, &surface_normal);
    return if emit { obj.emission + e } else { e }
        + obj.color.mul_element_wise(radiance(
                &config,
                &Ray {
                    origin: *intersect_point,
                    direction: random_direction,
                },
                depth + 1,
                rng,
                new_emit,
            ));
}

fn spec_radiance(
    config: &Config,
    ray: &Ray,
    depth: i32,
    rng: &mut ThreadRng,
    surface_normal: &Vec3d,
    intersect_point: &Vec3d,
    obj: &Sphere,
) -> Vec3d {
    obj.emission
        + obj.color.mul_element_wise(
            radiance(
                &config,
                &get_reflected_ray(&ray, &surface_normal, &intersect_point),
                depth + 1,
                rng,
                true,
            ))
}

fn refr_radiance(
    config: &Config,
    ray: &Ray,
    depth: i32,
    rng: &mut ThreadRng,
    surface_normal: &Vec3d,
    intersect_point: &Vec3d,
    oriented_surface_normal: &Vec3d,
    obj: &Sphere,
) -> Vec3d {
    let refl_ray = get_reflected_ray(&ray, &surface_normal, &intersect_point);
    let into = surface_normal.dot(*oriented_surface_normal) > 0.0;

    let air = 1.0;
    let glass = 1.5;
    let refraction = if into { air / glass } else { glass / air };

    let angle = ray.direction.dot(*oriented_surface_normal);
    let cos2t = 1.0 - refraction * refraction * (1.0 - angle * angle);
    if cos2t < 0.0 {
        // Total internal reflection so all light is reflected (internally)
        return obj.emission
            + obj.color.mul_element_wise(radiance(&config, &refl_ray, depth + 1, rng, true));
    }

    let refract_dir = (ray.direction * refraction
        - *surface_normal * (if into { 1.0 } else { -1.0 } * (angle * refraction + cos2t.sqrt())))
        .normalize();

    // Fresnel reflectance
    let a = glass - air;
    let b = glass + air;
    let normal_reflected = a * a / (b * b);
    let c = if into {
        1.0 + angle
    } else {
        1.0 - refract_dir.dot(*surface_normal)
    };

    let total_reflected = normal_reflected + (1.0 - normal_reflected) * c.powi(5);
    let total_refracted = 1.0 - total_reflected;
    let reflect_probability = 0.25 + 0.5 * total_reflected;

    // Weight results based on probability
    let reflect_weight = total_reflected / reflect_probability;
    let refract_weight = total_refracted / (1.0 - reflect_probability);

    return obj.emission
        + obj.color.mul_element_wise(
            // if depth is shallow, we want to sample everything
            if depth > 2 {
                if rng.next_f64() < reflect_probability {
                    // reflect ray
                    radiance(&config, &refl_ray, depth + 1, rng, true) * reflect_weight
                } else {
                    // refract ray
                    radiance(
                        &config,
                        &Ray {
                            origin: *intersect_point,
                            direction: refract_dir,
                        },
                        depth + 1,
                        rng,
                        true,
                    ) * refract_weight
                }
            } else {
                // do both
                radiance(&config, &refl_ray, depth + 1, rng, true) * total_reflected
                    + radiance(
                        &config,
                        &Ray {
                            origin: *intersect_point,
                            direction: refract_dir,
                        },
                        depth + 1,
                        rng,
                        true,
                    ) * total_refracted
            },
        );
}

pub fn radiance(config: &Config, ray: &Ray, depth: i32, rng: &mut ThreadRng, emit: bool) -> Vec3d {
    match config.scene.intersect(ray) {
        None => Vec3d::new(0.0, 0.0, 0.0),
        Some((closest_intersect_distance, id)) => {

            let obj: &Sphere = &config.scene.spheres[id]; // the hit object

            // Russian Roulette
            // Use maximum component (r,g,b) of the surface color
            let brightest_color = get_brightest_color(&obj.color); // brightest colour

            // Don't do Russian Roulette until after config.roulette_depth
            // More likely to end on a darker surface
            if brightest_color == 0.0
                || (depth > config.roulette_depth && rng.next_f64() >= brightest_color)
            {
                return if emit { obj.emission } else { Vec3d::new(0.0, 0.0, 0.0) };
            }

            let intersect_point = ray.origin + ray.direction * closest_intersect_distance; // point on sphere where intersects
            let surface_normal = (intersect_point - obj.position).normalize(); // surface oriented_surface_normal of intersection point


            match obj.reflection {
                ReflectType::DIFF => diff_radiance(
                    &config,
                    &ray,
                    depth,
                    rng,
                    emit,
                    &surface_normal,
                    &intersect_point,
                    &obj,
                ),
                ReflectType::SPEC => spec_radiance(
                    &config,
                    &ray,
                    depth,
                    rng,
                    &surface_normal,
                    &intersect_point,
                    &obj,
                ),
                ReflectType::REFR => {
                    // corrected oriented_surface_normal (ie internal or external intersection)
                    let oriented_surface_normal = if surface_normal.dot(ray.direction) < 0.0 {
                        // dot product negative if ray is internal
                        surface_normal
                    } else {
                        surface_normal * -1.0
                    };

                    refr_radiance(
                        &config,
                        &ray,
                        depth,
                        rng,
                        &surface_normal,
                        &intersect_point,
                        &oriented_surface_normal,
                        &obj,
                    )
                },
            }
        }
    }
}

#[cfg(all(feature = "unstable", test))]
mod bench {
    extern crate rand;
    extern crate test;

    use super::*;
    use self::test::Bencher;

    #[bench]
    fn bench_radiance(b: &mut Bencher) {
        let scene = Scene::new2();
        let mut rng = rand::thread_rng();
        let ray = Ray {
            origin: Vec3d::new( 50.0, 50.0, 100.0),
            direction: Vec3d::new( 0.0, -0.042612, -1.0).normalize(),
        };

        b.iter(|| radiance(&config, &ray, 0, &mut rng));
    }

    #[bench]
    fn bench_radiance_black(b: &mut Bencher) {
        let scene = Scene::new2();
        let mut rng = rand::thread_rng();
        let ray = Ray {
            origin: Vec3d::new( 50.0, 50.0, 260.0),
            direction: Vec3d::new( 0.0, -0.042612, -1.0).normalize(),
        };

        b.iter(|| radiance(&config, &ray, 0, &mut rng));
    }
}
