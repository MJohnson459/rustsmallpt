extern crate image;
extern crate rand;
extern crate rayon;

use std::io::{self, Write};
use image::ImageBuffer;
use rand::Rng;
use rayon::prelude::*;
use cgmath::InnerSpace;

use ray::Ray;
use ray::radiance;
use utility::*;
use vector_3d::{Vec3d, vec_clamp};
use config::Config;

pub struct Camera {
    ray: Ray,
    cx: Vec3d,
    cy: Vec3d,
}

impl Camera {
    pub fn new(config: &Config) -> Camera {
        let fov = 0.6;
        let position = Vec3d::new( 50.0, 50.0, 260.0);

        let ray: Ray = Ray {
            origin: position,
            direction: Vec3d::new(0.0, -0.042612, -1.0).normalize(),
        };
        let cx: Vec3d = Vec3d::new( (config.width as f64) * fov / (config.height as f64), 0.0, 0.0); // x direction increment
        let cy: Vec3d = (cx.cross(ray.direction)).normalize() * fov; // y direction increment

        Camera {
            ray: ray,
            cx: cx,
            cy: cy,
        }
    }

    fn update_pixel(&self, x: f64, y: f64, config: &Config, weight: f64) -> Vec3d {
        let width = config.width as f64;
        let height = config.height as f64;
        let mut rng = rand::thread_rng();

        let mut new_value: Vec3d = Vec3d::new( 0.0, 0.0, 0.0);

        if config.tent_filter {
            for sample_y in 0..2 {
                for sample_x in 0..2 {
                    let x_offset = calculate_offset(sample_x as f64);
                    let y_offset = calculate_offset(sample_y as f64);

                    let d = self.cx * ((x + x_offset) / width - 0.5)
                        + self.cy * ((y + y_offset) / height - 0.5)
                        + self.ray.direction;
                    let rad = radiance(
                        &config,
                        &Ray {
                            origin: self.ray.origin + d,
                            direction: d.normalize(),
                        },
                        0,
                        &mut rng,
                        true,
                    );

                    new_value += rad * 0.25; // 2x2 tent filter so weight is 0.25 each
                }
            }
        } else {
            let d: Vec3d =
                self.cx * (x / width - 0.5) + self.cy * (y / height - 0.5) + self.ray.direction;
            new_value = radiance(
                &config,
                &Ray {
                    origin: self.ray.origin + d,
                    direction: d.normalize(),
                },
                0,
                &mut rng,
                true,
            );
        }

        if config.explicit_light_sampling {
            vec_clamp(&new_value) * weight
        } else {
            new_value * weight
        }
    }

    fn single_sample(&self, prev_screen: &mut Vec<Vec3d>, config: &Config, weight: f64) {
        prev_screen
            .par_iter_mut()
            .enumerate()
            .for_each(|(i, value)| {
                let x = i % config.width;
                let y = i / config.width;
                *value = *value + self.update_pixel(x as f64, y as f64, config, weight);
            });
    }

    pub fn render_scene(&self, config: &Config) {
        let weight = 1.0 / config.samples as f64;
        let mut screen = vec![Vec3d::new(0.0, 0.0, 0.0); config.width * config.height];
        for sample in 0..config.samples {
            print!("Rendering {:.4}%\r", 100 * sample / config.samples);
            match io::stdout().flush() {
                Ok(_v) => {}
                Err(e) => println!("{}", e),
            }
            self.single_sample(&mut screen, &config, weight);
            if config.save_per_sample || sample == config.samples - 1 {
                let image = to_image(config.width, config.height, &screen);
                match image.save(config.path.as_path()) {
                    Ok(_v) => {}
                    Err(e) => println!("{}", e),
                }
            }
        }
    }
}

fn calculate_offset(bias: f64) -> f64 {
    let mut rng = rand::thread_rng();
    // Tent filter
    let dx: f64;
    let x = 2.0 * rng.next_f64();
    if x < 1.0 {
        dx = x.sqrt() - 1.0;
    } else {
        dx = 1.0 - (2.0 - x).sqrt();
    }

    (0.5 + dx + bias) / 2.0
}

fn to_image(
    width: usize,
    height: usize,
    screen: &Vec<Vec3d>,
) -> ImageBuffer<image::Rgb<u8>, Vec<u8>> {
    assert!(height * width == screen.len());

    ImageBuffer::from_fn(width as u32, height as u32, |x, y| {
        let i: usize = x as usize + (height - 1 - y as usize) * width;
        image::Rgb([to_u8(screen[i].x), to_u8(screen[i].y), to_u8(screen[i].z)])
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use scene::AvailableScenes;

    #[test]
    fn test_single_sample() {
        let width = 100;
        let height = 100;
        let weight = 1.0;
        let config = Config::new(width, height, 1, Some(AvailableScenes::Floating));
        let camera = Camera::new(&config);
        let mut prev_screen = vec![Vec3d::new(0.0, 0.0, 0.0); width * height];
        camera.single_sample(&mut prev_screen, &config, weight);
    }
}

#[cfg(all(feature = "unstable", test))]
mod bench {
    extern crate test;

    use super::*;
    use self::test::Bencher;
    use scene::AvailableScenes;

    #[bench]
    fn bench_single_sample(b: &mut Bencher) {
        let width = 100;
        let height = 100;
        let config = Config::new(width, height, 1, Some(AvailableScenes::Floating));
        let camera = Camera::new(&config);
        b.iter(|| {
            let mut prev_screen = vec![Vec3d::new(0.0, 0.0, 0.0); width * height];
            camera.single_sample(&mut prev_screen, &config, 1.0)
        });
    }

    #[bench]
    fn bench_update_pixel(b: &mut Bencher) {
        let width = 100;
        let height = 100;
        let camera = Camera::new(width, height);
        b.iter(|| camera.update_pixel(5.0, 5.0, width as f64, height as f64, 1.0));
    }

}
