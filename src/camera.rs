#![cfg_attr(feature = "unstable", feature(test))]

extern crate image;
extern crate rand;
extern crate rayon;

use std::io::{self, Write};
use image::ImageBuffer;
use rand::Rng;
use rayon::prelude::*;
use std::path::Path;
use na::{Vector3};

use ray::Ray;
use ray::radiance;
use scene::Scene;
use utility::*;

pub struct Camera {
    ray: Ray,
    cx: Vector3<f64>,
    cy: Vector3<f64>,
}

impl Camera {
    pub fn new(width: usize, height: usize) -> Camera {
        let fov = 0.4;
        let position = Vector3::<f64>::new(50.0, 50.0, 330.0);

        let ray: Ray = Ray{origin: position, direction: Vector3::<f64>::new(0.0, -0.042612, -1.0).normalize()};
        let cx = Vector3::<f64>::new((width as f64)*fov/(height as f64), 0.0, 0.0); // x direction increment
        let cy = (cx.cross(&ray.direction)).normalize() * fov;                    // y direction increment

        Camera {
            ray: ray,
            cx: cx,
            cy: cy,
        }
    }

    fn update_pixel(&self, prev_value: &Vector3<f64>, x: f64, y: f64, width: f64, height: f64, weight: f64, scene: &Scene) -> Vector3<f64> {
        let mut rng = rand::thread_rng();
        let tent_filter = true;
        let camera_ray_offset = 140.0;

        let mut new_value: Vector3<f64> = Vector3::<f64>::zeros();

        if tent_filter {
            for sample_y in 0..2 {
                for sample_x in 0..2 {
                    let x_offset = calculate_offset(sample_x as f64);
                    let y_offset = calculate_offset(sample_y as f64);

                    let mut d: Vector3<f64> = self.cx * ((x + x_offset) / width - 0.5) + self.cy * ((y + y_offset) / height - 0.5) + self.ray.direction;
                    let rad: Vector3<f64> = radiance(&scene, &Ray{origin: self.ray.origin + d * camera_ray_offset, direction: d.normalize()}, 0, &mut rng, 1.0);

                    new_value += rad * 0.25; // 2x2 tent filter so weight is 0.25 each
                }
            }
        } else {
            let mut d: Vector3<f64> = self.cx * (x / width - 0.5) + self.cy * (y / height - 0.5) + self.ray.direction;
            new_value = radiance(&scene, &Ray{origin: self.ray.origin + d * camera_ray_offset, direction: d.normalize()}, 0, &mut rng, 1.0);
        }


        *prev_value + (new_value * weight)
    }

    fn single_sample(&self, prev_screen: &Vec<Vector3<f64>> , width: usize, height: usize, weight: f64, scene: &Scene) -> Vec<Vector3<f64>> {
        let pixels = height * width;
        (0..pixels).into_par_iter().map(|i| {
            let x = i % width;
            let y = i / width;
            self.update_pixel(&prev_screen[i], x as f64, y as f64, width as f64, height as f64, weight, &scene)
        }).collect()
    }

    pub fn render_scene(&self, scene: &Scene, width: usize, height: usize, samples: usize, path: &Path) {
        let save_per_sample = true;
        let weight = 1.0 / samples as f64;
        let mut prev_screen = vec![Vector3::<f64>::zeros(); width * height];
        for sample in 0..samples {
            print!("Rendering {:.4}%\r", 100 * sample / samples);
            io::stdout().flush().unwrap();
            let new_screen = self.single_sample(&prev_screen, width, height, weight, &scene);
            if save_per_sample || sample == samples - 1 {
                let image = to_image(width, height, &new_screen);
                image.save(&path).expect("Unable to save image at path");
            }
            prev_screen = new_screen;
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

fn to_image(width: usize, height: usize, screen: &Vec<Vector3<f64>>) -> ImageBuffer<image::Rgb<u8>, Vec<u8>> {
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
        let prev_screen = vec!(Vector3<f64>::zeros(); width * height);
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
