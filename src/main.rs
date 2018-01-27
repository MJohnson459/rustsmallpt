#![cfg_attr(feature = "unstable", feature(test))]

#[macro_use]
extern crate serde_derive;
extern crate docopt;
extern crate time;
extern crate image;
extern crate rand;
extern crate rayon;

use docopt::Docopt;
use std::path::Path;

mod camera;
mod material;
mod ray;
mod scene;
mod sphere;
mod utility;
mod vector_3d;

use camera::Camera;
use scene::AvailableScenes;
use scene::Scene;
use utility::*;


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

fn print_estimate(width: usize, height: usize, samples: usize) {
    println!("width: {}, height: {}, samples: {}", width, height, samples);
    println!("Number of threads: {}", rayon::current_num_threads());

    let time_per_spp: f64 = 1.8612e-6;
    let est_time: f64 = time_per_spp*(samples*width*height) as f64;
    println!("Estimated time [RELEASE]: {}", format_time(est_time));

    let time_per_spp: f64 = 2.0857e-6;
    let est_time: f64 = time_per_spp*(samples*width*height) as f64;
    println!("Estimated time [RELEASE][SAVE]: {}", format_time(est_time));
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
                AvailableScenes::New1 => scene = Scene::new(),
                AvailableScenes::New2 => scene = Scene::new2(),
            }
        },
        None => {
            println!("Using default scene New2");
            scene = Scene::new2();
        }
    }

    let time_start = time::precise_time_s();


    let image_name = format!("{}_{}_{}_{}.png", scene.name, width, height, samples);
    let path = Path::new(&image_name);

    let camera = Camera::new(width, height);
    camera.render_scene(&scene, width, height, samples, path);

    let time_taken = time::precise_time_s() - time_start;

    println!("Finished rendering. Time taken: {}", format_time(time_taken));
    println!("time_per_spp*1e6: {}", (time_taken as f64/(width*height*samples) as f64)*1e6);
    println!("samples per second: {}", (width*height*samples) as f64/time_taken);
}

