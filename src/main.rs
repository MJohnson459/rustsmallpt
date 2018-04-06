#![cfg_attr(feature = "unstable", feature(test))]

extern crate docopt;
extern crate image;
extern crate rand;
extern crate rayon;
#[macro_use] extern crate serde_derive;
extern crate time;
extern crate cgmath;
extern crate num_traits;

use docopt::Docopt;

mod camera;
mod config;
mod material;
mod ray;
mod scene;
mod sphere;
mod utility;
mod vector_3d;

use camera::Camera;
use config::Config;
use scene::AvailableScenes;
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
                    Valid values: Default, Floating, Lightbulb
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

    let time_per_spp: f64 = 1.7612e-6;
    let est_time: f64 = time_per_spp * (samples * width * height) as f64;
    println!("Estimated time [RELEASE]: {}", format_time(est_time));

    let time_per_spp: f64 = 1.9857e-6;
    let est_time: f64 = time_per_spp * (samples * width * height) as f64;
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

    let time_start = time::precise_time_s();

    let config = Config::new(width, height, samples, args.flag_scene);

    let camera = Camera::new(&config);
    camera.render_scene(&config);

    let time_taken = time::precise_time_s() - time_start;

    println!(
        "Finished rendering. Time taken: {}",
        format_time(time_taken)
    );
    println!(
        "time_per_spp*1e6: {}",
        (time_taken as f64 / (width * height * samples) as f64) * 1e6
    );
    println!(
        "samples per second: {}",
        (width * height * samples) as f64 / time_taken
    );
}
