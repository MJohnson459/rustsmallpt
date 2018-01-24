#![cfg_attr(feature = "unstable", feature(test))]

extern crate image;
extern crate rand;
extern crate rayon;

use image::ImageBuffer;
use rand::Rng;
use rayon::prelude::*;

use ray::Ray;
use ray::radiance;
use scene::Scene;
use utility::*;
use vector_3d::Vec3d;

pub struct Camera {
    ray: Ray,
    cx: Vec3d,
    cy: Vec3d
}

impl Camera {
    pub fn new(width: usize, height: usize, fov: f64) -> Camera {

        let ray: Ray = Ray{origin: Vec3d{x:50.0,y:50.0,z:295.6}, direction: Vec3d{x:0.0,y:-0.042612, z:-1.0}.normalise()};
        let cx: Vec3d = Vec3d{x:(width as f64)*fov/(height as f64),y:0.0,z:0.0}; // x direction increment
        let cy: Vec3d = (cx % ray.direction).normalise()*fov;                    // y direction increment

        Camera {
            ray: ray,
            cx: cx,
            cy: cy,
        }
    }

    pub fn render_scene(&self, scene: &Scene, width: usize, height: usize, samples: usize) -> ImageBuffer<image::Rgb<u8>, Vec<u8>> {
        let screen = self.single_sample(width, height, samples, &scene);
        to_image(&screen)
    }

    fn update_pixel(&self, x: f64, y: f64, width: f64, height: f64, samples: usize, scene: &Scene) -> Vec3d {
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
                sum = sum + v*0.25;
            }
        }
        return sum;
    }

    fn single_sample(&self, width: usize, height: usize, samples: usize, scene: &Scene) -> Vec<Vec<Vec3d>> {
        // let img = RgbImage::new(width as u32, height as u32);
        (0..height).map(|y| {
            print!("Rendering ({} spp) {:.4}%...\r", samples * 4, 100.0 * y as f64 / height as f64);
            (0..width).into_par_iter().map(|x| {
                self.update_pixel(x as f64, y as f64, width as f64, height as f64, samples, &scene)
            }).collect()
        }).collect()
    }

}



fn to_image(screen: &Vec<Vec<Vec3d>>) -> ImageBuffer<image::Rgb<u8>, Vec<u8>> {
    let height = screen.len();
    let width = screen[0].len();

    ImageBuffer::from_fn(width as u32, height as u32, |x, y| {
        let w = x as usize;
        let h = height - 1 - y as usize;
        image::Rgb([to_u8(screen[h][w].x), to_u8(screen[h][w].y), to_u8(screen[h][w].z)])
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
