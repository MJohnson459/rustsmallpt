#![cfg_attr(feature = "unstable", feature(test))]

#[macro_use]
extern crate serde_derive;

extern crate docopt;
extern crate image;
extern crate rand;
extern crate rayon;
extern crate time;

use docopt::Docopt;
use rand::Rng;
use rayon::prelude::*;
use std::path::Path;

mod camera;
mod material;
mod ray;
mod scene;
mod sphere;
mod utility;
mod vector_3d;

use ray::Ray;
use scene::AvailableScenes;
use scene::Scene;
use utility::*;
use vector_3d::Vec3d;


// Write the Docopt usage string.
static USAGE: &'static str = "
Rust (not-so-)small Path Tracer

Usage:
    smallpt <width> <height> <samples> [--scene SCENE]
    smallpt (--help | --version)

Options:
    -h --help       Show this screen
    --scene SCENE   The scene to render
                    Valid values: New1, New2
";

#[derive(Debug, Deserialize)]
struct Args {
    arg_width: usize,
    arg_height: usize,
    arg_samples: usize,
    flag_scene: Option<AvailableScenes>,
    flag_help: bool,
    flag_version: bool,
}

fn update_pixel(x: f64, y: f64, width: f64, height: f64, samples: usize, scene: &Scene) -> Vec3d {
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

                let mut d: Vec3d = scene.cam.cx*((((sample_x as f64)+0.5 + dx)/2.0 + x) / width - 0.5) +
                                scene.cam.cy*((((sample_y as f64)+0.5 + dy)/2.0 + y) / height - 0.5) + scene.cam.ray.direction;
                let rad: Vec3d = scene.radiance(&Ray{origin: scene.cam.ray.origin + d*140.0, direction: d.normalise()}, 0, &mut rng, 1.0);
                r = r + rad*(1.0/samples as f64);
            }
            let v: Vec3d = Vec3d{x: clamp(r.x), y: clamp(r.y), z: clamp(r.z)};
            //println!("v.x: {}, v.y: {}, v.z: {}", v.x, v.y, v.z);
            sum = sum + v*0.25;
        }
    }
    return sum;
}

pub fn single_sample(width: usize, height: usize, samples: usize, scene: &Scene) -> Vec<Vec<Vec3d>> {
    // let img = RgbImage::new(width as u32, height as u32);
    (0..height).map(|y| {
        print!("Rendering ({} spp) {:.4}%...\r", samples * 4, 100.0 * y as f64 / height as f64);
        (0..width).into_par_iter().map(|x| {
            update_pixel(x as f64, y as f64, width as f64, height as f64, samples, &scene)
        }).collect()
    }).collect()
}

fn print_estimate(width: usize, height: usize, samples: usize) {
    println!("width: {}, height: {}, samples: {}", width, height, samples);
    println!("Number of threads: {}", rayon::current_num_threads());

    let time_per_spp: f64 = 3.659458e-6;
    let est_time: f64 = 4.0*time_per_spp*(samples*width*height) as f64;

    println!("Estimated time [ DEBUG ]: {}", format_time(est_time));

    let time_per_spp: f64 = 0.35e-6;
    let est_time: f64 = 4.0*time_per_spp*(samples*width*height) as f64;

    println!("Estimated time [RELEASE]: {}", format_time(est_time));
}


fn main() {
    // Parse argv and exit the program with an error message if it fails.
    let args: Args = Docopt::new(USAGE)
        .and_then(|dopt| dopt.deserialize())
        .unwrap_or_else(|e| e.exit());

    let width: usize = args.arg_width;
    let height: usize = args.arg_height;
    let samples: usize = args.arg_samples;

    print_estimate(width, height, samples);

    let scene: Scene;
    match args.flag_scene {
        Some(x) => {
            match x {
                AvailableScenes::New1 => scene = Scene::new(width, height),
                AvailableScenes::New2 => scene = Scene::new2(width, height),
            }
        },
        None => {
            println!("Using default scene New2");
            scene = Scene::new2(width, height);
        }
    }

    let time_start = time::precise_time_s();

    let screen = single_sample(width, height, samples, &scene);
    let image = to_image(&screen);

    let time_taken = time::precise_time_s() - time_start;

    println!("Finished rendering. Time taken: {}", format_time(time_taken));
    println!("DEBUG time_per_spp*1e6: {}", (time_taken as f64/(width*height*4*samples) as f64)*1e6);

    let image_name = format!("{}_{}_{}_{}.png", scene.name, width, height, samples*4);
    image.save(&Path::new(&image_name)).unwrap();
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
