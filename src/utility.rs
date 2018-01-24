#![cfg_attr(feature = "unstable", feature(test))]

extern crate image;
extern crate time;

use image::ImageBuffer;

use vector_3d::Vec3d;

pub fn float_eq(a: f64, b:f64) -> bool {
    (a - b).abs() < 0.0001
}

pub fn clamp(x: f64) -> f64 {
    if x < 0.0 {
        0.0
    } else if x > 1.0 {
        1.0
    } else {
        x
    }
}

pub fn to_int(x: f64) -> i32 {
    (clamp(x).powf(1.0/2.2)*255.0+0.5) as i32
}

pub fn to_u8(x: f64) -> u8 {
    (clamp(x).powf(1.0/2.2)*255.0+0.5) as u8
}

pub fn format_time(seconds: f64) -> String {
    let hours: f64 = (seconds/3600.0).floor();
    let mins: f64 = ((seconds - hours*3600.0)/60.0).floor();
    let secs: f64 = seconds - mins*60.0;
    format!("{:02.0}:{:02.0}:{:02.0}", hours, mins, secs)

}

pub fn to_image(screen: &Vec<Vec<Vec3d>>) -> ImageBuffer<image::Rgb<u8>, Vec<u8>> {
    let height = screen.len();
    let width = screen[0].len();

    ImageBuffer::from_fn(width as u32, height as u32, |x, y| {
        let w = x as usize;
        let h = height - 1 - y as usize;
        image::Rgb([to_u8(screen[h][w].x), to_u8(screen[h][w].y), to_u8(screen[h][w].z)])
    })
}