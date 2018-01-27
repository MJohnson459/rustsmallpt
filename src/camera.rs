#![cfg_attr(feature = "unstable", feature(test))]

extern crate image;
extern crate rand;
extern crate rayon;

use std::io::{self, Write};
use image::ImageBuffer;
use rand::Rng;
use rayon::prelude::*;
use std::path::Path;

use ray::Ray;
use ray::radiance;
use scene::Scene;
use utility::*;
use vector_3d::Vec3d;

pub struct Camera {
    ray: Ray,
    cx: Vec3d,
    cy: Vec3d,
}

impl Camera {
    pub fn new(width: usize, height: usize) -> Camera {
        let fov = 0.4;
        let position = Vec3d{x: 50.0, y: 50.0, z: 330.0};

        let ray: Ray = Ray{origin: position, direction: Vec3d{x: 0.0, y: -0.042612, z: -1.0}.normalise()};
        let cx: Vec3d = Vec3d{x:(width as f64)*fov/(height as f64),y:0.0,z:0.0}; // x direction increment
        let cy: Vec3d = (cx % ray.direction).normalise()*fov;                    // y direction increment

        Camera {
            ray: ray,
            cx: cx,
            cy: cy,
        }
    }

    fn update_pixel(&self, x: f64, y: f64, width: f64, height: f64, weight: f64, scene: &Scene) -> Vec3d {
        let mut rng = rand::thread_rng();
        let tent_filter = true;
        let camera_ray_offset = 140.0;

        let mut new_value: Vec3d = Vec3d{x:0.0, y: 0.0, z: 0.0};

        if tent_filter {
            for sample_y in 0..2 {
                for sample_x in 0..2 {
                    let x_offset = calculate_offset(sample_x as f64);
                    let y_offset = calculate_offset(sample_y as f64);

                    let mut d: Vec3d = self.cx * ((x + x_offset) / width - 0.5) + self.cy * ((y + y_offset) / height - 0.5) + self.ray.direction;
                    let rad: Vec3d = radiance(&scene, &Ray{origin: self.ray.origin + d * camera_ray_offset, direction: d.normalise()}, 0, &mut rng, 1.0);

                    new_value += rad * 0.25; // 2x2 tent filter so weight is 0.25 each
                }
            }
        } else {
            let mut d: Vec3d = self.cx * (x / width - 0.5) + self.cy * (y / height - 0.5) + self.ray.direction;
            new_value = radiance(&scene, &Ray{origin: self.ray.origin + d * camera_ray_offset, direction: d.normalise()}, 0, &mut rng, 1.0);
        }

        new_value * weight
    }

    fn single_sample(&self, prev_screen: &mut Vec<Vec3d> , width: usize, height: usize, weight: f64, scene: &Scene) {
        prev_screen.par_iter_mut().enumerate().for_each(|(i, value)| {
            let x = i % width;
            let y = i / width;
            *value = *value + self.update_pixel(x as f64, y as f64, width as f64, height as f64, weight, &scene);
        });
    }

    pub fn render_scene(&self, scene: &Scene, width: usize, height: usize, samples: usize, path: &Path) {
        let save_per_sample = true;
        let weight = 1.0 / samples as f64;
        let mut screen = vec![Vec3d::default(); width * height];
        for sample in 0..samples {
            print!("Rendering {:.4}%\r", 100 * sample / samples);
            io::stdout().flush().unwrap();
            self.single_sample(&mut screen, width, height, weight, &scene);
            if save_per_sample || sample == samples - 1 {
                let image = to_image(width, height, &screen);
                image.save(&path).expect("Unable to save image at path");
            }
        }
    }

}

fn calculate_offset(bias: f64) -> f64 {
    let mut rng = rand::thread_rng();
    // Tent filter
    let dx: f64;
    let x = 2.0 * rng.next_f64(); //erand48(xi);
    if x < 1.0 {
        dx = x.sqrt() - 1.0;
    } else {
        dx = 1.0 - (2.0 - x).sqrt();
    }

    (0.5 + dx + bias) / 2.0
}

fn to_image(width: usize, height: usize, screen: &Vec<Vec3d>) -> ImageBuffer<image::Rgb<u8>, Vec<u8>> {
    assert!(height * width == screen.len());

    ImageBuffer::from_fn(width as u32, height as u32, |x, y| {
        let i: usize = x as usize + (height - 1 - y as usize) * width;
        image::Rgb([to_u8(screen[i].x), to_u8(screen[i].y), to_u8(screen[i].z)])
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_sample() {
        let width = 100;
        let height = 100;
        let weight = 1.0;
        let scene = Scene::new2();
        let camera = Camera::new(width, height);
        let prev_screen = vec!(Vec3d::default(); width * height);
        let result = camera.single_sample(&prev_screen, width, height, weight, &scene);
        assert_eq!(result.len(), prev_screen.len());
    }
}

#[cfg(all(feature = "unstable", test))]
mod bench {
    extern crate test;

    use super::*;
    use self::test::Bencher;

    #[bench]
    fn bench_single_sample(b: &mut Bencher) {
        let width = 100;
        let height = 100;
        let samples = 2;
        let scene = Scene::new2(width, height);
        b.iter(|| single_sample(width, height, samples, &scene));
    }

    #[bench]
    fn bench_update_pixel(b: &mut Bencher) {
        let width = 100;
        let height = 100;
        let samples = 2;
        let scene = Scene::new2(width, height);
        b.iter(|| update_pixel(5.0, 5.0, width as f64, height as f64, samples, &scene));
    }

}
