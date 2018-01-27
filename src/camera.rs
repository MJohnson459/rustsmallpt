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
    width: usize,
    height: usize,
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
            width: width,
            height: height,
            ray: ray,
            cx: cx,
            cy: cy,
        }
    }

//    pub fn render_scene(&self, scene: &Scene, width: usize, height: usize, samples: usize) -> ImageBuffer<image::Rgb<u8>, Vec<u8>> {
//        let screen = self.single_sample(width, height, samples, &scene);
//        to_image(width, height, &screen)
//    }

    fn update_pixel_tent(&self, x: f64, y: f64, width: f64, height: f64, samples: usize, scene: &Scene) -> Vec3d {
        let mut rng = rand::thread_rng();
        let mut sum = Vec3d{x:0.0,y:0.0,z:0.0};
        for sample_y in 0..2 {
            for sample_x in 0..2 {
                let mut r: Vec3d = Vec3d{x:0.0, y: 0.0, z: 0.0};
                for _ in 0..samples {
                    let rand_x = 2.0*rng.next_f64(); //erand48(xi);
                    let rand_y = 2.0*rng.next_f64(); //erand48(xi);

                    let dx: f64;
                    let dy: f64;

                    // Tent filter
                    // rand_x and rand_y determine location within pixel
                    if rand_x < 1.0 {
                        dx = rand_x.sqrt() - 1.0;
                    } else {
                        dx = 1.0 - (2.0 - rand_x).sqrt();
                    }

                    if rand_y < 1.0 {
                        dy = rand_y.sqrt() - 1.0;
                    } else {
                        dy = 1.0 - (2.0 - rand_y).sqrt();
                    }

                    let mut d: Vec3d = self.cx*((((sample_x as f64)+0.5 + dx)/2.0 + x) / width - 0.5) +
                                    self.cy*((((sample_y as f64)+0.5 + dy)/2.0 + y) / height - 0.5) + self.ray.direction;
                    let rad: Vec3d = radiance(&scene, &Ray{origin: self.ray.origin + d*140.0, direction: d.normalise()}, 0, &mut rng, 1.0);
                    r = r + rad*(1.0/samples as f64);
                }
                let v: Vec3d = Vec3d{x: clamp(r.x), y: clamp(r.y), z: clamp(r.z)};
                //println!("v.x: {}, v.y: {}, v.z: {}", v.x, v.y, v.z);
                sum += v*0.25;
            }
        }
        return sum;
    }

    fn update_pixel(&self, prev_value: &Vec3d, x: f64, y: f64, width: f64, height: f64, weight: f64, scene: &Scene) -> Vec3d {
        let mut rng = rand::thread_rng();
        let mut d: Vec3d = self.cx * (x / width - 0.5) + self.cy * (y / height - 0.5) + self.ray.direction;
        let rad: Vec3d = radiance(&scene, &Ray{origin: self.ray.origin + d*140.0, direction: d.normalise()}, 0, &mut rng, 1.0);

        *prev_value + (rad * weight)
    }

    fn single_sample(&self, prev_screen: &Vec<Vec3d> , width: usize, height: usize, weight: f64, scene: &Scene) -> Vec<Vec3d> {
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
        let mut prev_screen = vec![Vec3d::default(); width * height];
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
        let samples = 1;
        let scene = Scene::new2(width, height);
        let result = single_sample(width, height, samples, &scene);
        assert_eq!(result.len(), height);
        assert_eq!(result[0].len(), width);
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
