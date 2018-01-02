#![cfg_attr(feature = "unstable", feature(test))]

extern crate rand;
extern crate threadpool;
extern crate num_cpus;
extern crate time;
extern crate docopt;
extern crate image;

use std::env;
use std::f64::consts::PI;
use std::io::Write;
use std::io;
use std::sync::Arc;
use threadpool::ThreadPool;
use std::sync::mpsc::channel;
use docopt::Docopt;
use rand::Rng;
use rand::ThreadRng;

mod vector_3d;

use vector_3d::Vec3d;

use std::path::Path;

// Write the Docopt usage string.
static USAGE: &'static str = "
Usage:
    smallpt <width> <height> <samples>
";

pub fn float_eq(a: f64, b:f64) -> bool {
    (a - b).abs() < 0.0001
}

// ----------- Ray --------------

#[derive(Copy, Clone, Debug)]
pub struct Ray {
    origin: Vec3d,
    direction: Vec3d,
}

// ----------
#[derive(Copy, Clone, Debug)]
pub enum ReflectType {
    DIFF,
    SPEC,
    REFR
}

// ------
#[derive(Copy, Clone, Debug)]
pub struct Sphere {
    radius: f64,
    position: Vec3d,
    emission: Vec3d,
    color: Vec3d,
    reflection: ReflectType
}

impl Sphere {
    fn intersect(&self, ray: &Ray) -> f64 {
        let op: Vec3d = self.position - ray.origin;
        let eps: f64 = 1e-4;
        let b: f64 = op.dot(&ray.direction);
        let mut det: f64 = b*b - op.dot(&op) + self.radius*self.radius;
        if det < 0.0 {
            //println!("det < 0.0 : {}", det);
            return 0.0
        } else {
            //println!("det >= 0.0 : {}", det);
            det = det.sqrt();
        }

        //println!("\nself: {:?}",self);
        //println!("ray: {:?}",ray);
        //println!("op: {:?}",op);
        //println!("eps: {:?}",eps);
        //println!("b: {:?}",b);
        //println!("det: {:?}",det);


        let t1 = b - det;
        let t2 = b + det;
        if t1 > eps {
            return t1;
        } else if t2 > eps {
            return t2;
        } else {
            return 0.0;
        }
    }
}

// ---------------------

pub fn clamp(x: f64) -> f64 {
    if x < 0.0 {
        return 0.0;
    } else     if x > 1.0 {
        return 1.0;
    }
    x
}

pub fn to_int(x: f64) -> i32 {
    (clamp(x).powf(1.0/2.2)*255.0+0.5) as i32
}

pub fn to_u8(x: f64) -> u8 {
    (clamp(x).powf(1.0/2.2)*255.0+0.5) as u8
}

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

        return Camera {
            ray: ray,
            cx: cx,
            cy: cy
        };
    }
}

pub struct Scene {
    pub spheres: Vec<Sphere>,
    pub cam: Camera
}

impl Scene {
    pub fn new(width: usize, height: usize) -> Scene {
        let mut spheres = Vec::new();
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:1e5+1.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.25,z:0.25}, reflection: ReflectType::DIFF}); // left
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:-1e5+99.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.25,y:0.25,z:0.75}, reflection: ReflectType::DIFF}); // right
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:40.8,z:1e5}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // back
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:40.8,z:-1e5+170.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // front
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:1e5,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // bottom
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:-1e5+81.6,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // Top
        spheres.push(Sphere {radius:16.5, position: Vec3d{x:27.0,y:16.5,z:47.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.999, reflection: ReflectType::SPEC}); // mirror
        spheres.push(Sphere {radius:16.5, position: Vec3d{x:73.0,y:16.5,z:78.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.999, reflection: ReflectType::REFR}); // glass
        spheres.push(Sphere {radius:600.0, position: Vec3d{x:50.0,y:681.6-0.27,z:81.6}, emission: Vec3d{x:12.0,y:12.0,z:12.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light

        Scene {
            spheres: spheres,
            cam: Camera::new(width, height, 0.5135)
        }
    }

    pub fn new2(width: usize, height: usize) -> Scene {
        let mut spheres = Vec::new();
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:1e5+1.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.8,y:0.1,z:0.1}, reflection: ReflectType::DIFF}); // left
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:-1e5+99.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.1,y:0.3,z:0.70}, reflection: ReflectType::DIFF}); // right
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:40.8,z:1e5}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // back
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:40.8,z:-1e5+170.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // front
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:1e5,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // bottom
        spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:-1e5+81.6,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // Top
        spheres.push(Sphere {radius:16.5, position: Vec3d{x:27.0,y:16.5,z:47.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.95, reflection: ReflectType::SPEC}); // mirror
        spheres.push(Sphere {radius:16.5, position: Vec3d{x:73.0,y:30.0,z:78.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.95, reflection: ReflectType::REFR}); // glass
        //spheres.push(Sphere {radius:600.0, position: Vec3d{x:50.0,y:681.6-0.27,z:81.6}, emission: Vec3d{x:4.0,y:4.0,z:4.0}*10000.0, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light
        spheres.push(Sphere {radius:600.0, position: Vec3d{x:50.0,y:681.6-0.27,z:81.6}, emission: Vec3d{x:12.0,y:12.0,z:12.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light
        //spheres.push(Sphere {radius:1.5, position: Vec3d{x:50.0,y:81.6-16.5,z:81.6}, emission: Vec3d{x:4.0,y:4.0,z:4.0}*100.0, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light
        Scene {
            spheres: spheres,
            cam: Camera::new(width, height, 0.5135)
        }
    }

    pub fn intersect(&self, ray: &Ray) -> (bool, f64, usize) {
        let n = self.spheres.len();
        let mut d: f64;
        let mut t: f64 = 1e20;
        let mut id: usize = 0;
        let inf = 1e20;

        for i in (0..n).rev() {
            //println!("i: {}",i);
            d = self.spheres[i].intersect(ray);
            if d != 0.0 && d < t {
                //println!("ray intersected sphere {}: {:?}", i, ray);
                t = d;
                id = i;
            }
        }
        //println!("intersects - t: {}, inf: {}", t, inf);
        (t < inf, t, id)
    }

    pub fn radiance(&self, ray: &Ray, mut depth: i32, rng: &mut ThreadRng, emission: f64) -> Vec3d {

        //println!("radiance - depth: {},  ray: {:?}", depth, ray);
        let (intersects, t, id) = self.intersect(ray);
        if !intersects {  // if miss, return black
            //println!("No Intersect!");
            return Vec3d{x:0.0, y:0.0, z:0.0};
        }

        let obj: Sphere = self.spheres[id].clone();  // the hit object

        let intersect_point: Vec3d = ray.origin + ray.direction*t; // point on sphere where intersects
        let surface_normal: Vec3d = (&intersect_point - &obj.position).normalise(); // surface light_normal of intersection point
        let light_normal: Vec3d; // corrected light_normal (ie internal or external intersection)
        if surface_normal.dot(&ray.direction) < 0.0 { // dot product negative if ray is internal
            light_normal = surface_normal;
        } else {
            light_normal = surface_normal*-1.0;
        }

        //println!("ray: {:?}",ray);
        //println!("t: {:?}",t);
        //println!("obj: {:?}",obj);
        //println!("x: {:?}",x);
        //println!("surface_normal: {:?}",surface_normal);
        //println!("light_normal: {:?}",light_normal);

        let mut f: Vec3d = obj.color.clone();

        let p: f64; // brightest colour
        if f.x > f.y && f.x>f.z {
            p = f.x;
        } else if f.y > f.z {
            p = f.y;
        } else {
            p = f.z;
        }

        depth = depth + 1;
        if depth > 5 || p == 0.0 {
            if rng.next_f64() < p {
                f = f*(1.0/p); // normalise colour [0,1]
            } else {
                // This might be at the wrong if statement
                return obj.emission*(emission as f64);
            }
        }

        match obj.reflection {
            ReflectType::DIFF => {
                let random_angle: f64 = 2.0*PI*rng.next_f64();
                let random: f64 = rng.next_f64();
                let r2s = random.sqrt();

                let w: Vec3d = light_normal.clone();
                let u: Vec3d;
                if w.x.abs() > 0.1 {
                    u = (Vec3d{x:0.0, y:1.0,z:0.0}.cross(w)).normalise();
                } else {
                    u = (Vec3d{x:1.0, y:0.0,z:0.0}.cross(w)).normalise();
                }
                assert!((u.length() - 1.0).abs() < 0.001);

                let v: Vec3d = w.cross(u);

                let random_direction: Vec3d = (u*random_angle.cos()*r2s + v*random_angle.sin()*r2s + w*(1.0-random).sqrt()).normalise();
                assert!((random_direction.length() - 1.0).abs() < 0.001);


                let e: Vec3d = Vec3d::zeros();
/*
                for i in 0..self.spheres.len() {
                    let sphere = &self.spheres[i];
                    if sphere.emission.x <= 0.0 &&
                        sphere.emission.y <= 0.0 &&
                        sphere.emission.z <= 0.0 { continue; };

                    let sw: Vec3d = sphere.position - intersect_point;
                    let su: Vec3d;
                    if sw.x.abs() > 0.1 {
                        su = (Vec3d::new(0.0,1.0,0.0)%sw).normalise();
                    } else {
                        su = (Vec3d::new(1.0,0.0,0.0)%sw).normalise();
                    }
                    let sv: Vec3d = sw%su;

                    let r = 40.0; //sphere.radius;

                    let cos_a_max: f64 = (1.0-r*r/(intersect_point-sphere.position).dot(&(intersect_point-sphere.position))).sqrt();
                    let eps1: f64 = rng.next_f64();
                    let eps2: f64 = rng.next_f64();
                    let cos_a = 1.0-eps1+eps1*cos_a_max;
                    let sin_a = (1.0-cos_a*cos_a).sqrt();
                    let phi = 2.0*PI*eps2;
                    let l: Vec3d = (su*phi.cos()*sin_a + sv*phi.sin()*sin_a + sw*cos_a).normalise();
                    let (intersects, _, id) = self.intersect(&Ray{origin: intersect_point, direction: l});
                    if intersects && id == i {
                        let omega: f64 = 2.0*PI*(1.0-cos_a_max);
                        e = e + f.mult(sphere.emission*l.dot(&light_normal)*omega)*FRAC_1_PI;
                    }
                }*/

                let fmu = f.mult(self.radiance(&Ray{origin: intersect_point, direction: random_direction}, depth, rng, 1.0));
                return obj.emission*(emission as f64) + e + fmu;
            },
            ReflectType::SPEC => {
                //println!("obj.emission SPEC, depth: {}", depth);
                let d: Vec3d = ray.direction-(surface_normal*2.0*surface_normal.dot(&ray.direction));
                return obj.emission + f.mult(self.radiance(&Ray{origin: intersect_point, direction: d}, depth, rng, 1.0));
            },
                _ => {}

        }

        let refl_ray: Ray = Ray{origin: intersect_point, direction: ray.direction-surface_normal*2.0*surface_normal.dot(&ray.direction)}; // Ideal dielectric REFRACTION
        let into: bool = surface_normal.dot(&light_normal) > 0.0;

        let nc: f64 = 1.0;
        let nt: f64 = 1.5;
        let nnt: f64;
        if into { // Ray from outside going in?
            nnt = nc/nt;
        } else {
            nnt = nt/nc;
        }
        let ddn: f64 =ray.direction.dot(&light_normal);
        let cos2t: f64 = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
        if cos2t < 0.0 { // Total internal reflection
            //println!("obj.emission cos2t < 0.0");
            return obj.emission + f.mult(self.radiance(&refl_ray, depth, rng, 1.0));
        }

        let tdir: Vec3d;

        if into {
            tdir = (ray.direction*nnt - surface_normal*(1.0*(ddn*nnt+cos2t.sqrt()))).normalise();
        } else {
            tdir = (ray.direction*nnt - surface_normal*(-1.0*(ddn*nnt+cos2t.sqrt()))).normalise();
        }

        let a: f64 = nt- nc;
        let b: f64 = nt + nc;
        let r0: f64 = a*a/(b*b);
        let c: f64;

        if into {
            c = 1.0 + ddn;
        } else {
            c = 1.0 - tdir.dot(&surface_normal);
        }

        let re: f64 = r0+(1.0-r0)*c.powi(5);
        let tr: f64 = 1.0 - re;
        let p: f64 = 0.25+0.5*re;
        let rp: f64 = re/p;
        let tp: f64 = tr/(1.0-p);

        if depth > 2 {
            if rng.next_f64() < p {
                //println!("rand::random::<f64>() < p - depth: {}", depth);
                return obj.emission + f.mult(self.radiance(&refl_ray, depth, rng, 1.0)*rp);
            } else {
                //println!("rand::random::<f64>() >= p - depth: {}", depth);
                return obj.emission + f.mult(self.radiance(&Ray{origin: intersect_point, direction: tdir}, depth, rng, 1.0)*tp);
            }
        } else {
            //println!("reflect ray: {}", depth);
            return obj.emission + f.mult(self.radiance(&refl_ray, depth, rng, 1.0)*re+self.radiance(&Ray{origin: intersect_point, direction: tdir}, depth, rng, 1.0)*tr);
        }
    }

}

struct Time {
    seconds: f64
}

impl Time {
    pub fn new(seconds: f64) -> Time {
        Time { seconds: seconds}
    }

    pub fn get_time(self) -> String {
        let hours: f64 = (self.seconds/3600.0).floor();
        let mins: f64 = ((self.seconds - hours*3600.0)/60.0).floor();
        let secs: f64 = self.seconds - mins*60.0;
        format!("{:02.0}:{:02.0}:{:02.0}", hours, mins, secs)

    }
}

fn update_pixel(x: f64, y: f64, width: f64, height: f64, samples: usize, scene: &Arc<Scene>) -> Vec3d {
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

pub fn single_row(y: usize, width: usize, height: usize, samples: usize, scene: Arc<Scene>) -> Vec<Vec3d> {
    let mut line = Vec::with_capacity(width);

    for x in 0..width {
        line.push(update_pixel(x as f64, y as f64, width as f64, height as f64, samples, &scene));
    }
    return line;
}


fn main() {

    // Prints each argument on a separate line
    for argument in env::args() {
        println!("{}", argument);
    }

    // Parse argv and exit the program with an error message if it fails.
    let args = Docopt::new(USAGE)
    .and_then(|dopt| dopt.parse())
    .unwrap_or_else(|e| e.exit());

    let width: usize = args.get_str("<width>").parse().unwrap_or(1024);
    let height: usize = args.get_str("<height>").parse().unwrap_or(768);
    let samples: usize = args.get_str("<samples>").parse().unwrap_or(100);


    println!("width: {}, height: {}, samples: {}", width, height, samples);

    let num_threads = num_cpus::get();

    let time_per_spp: f64 = 3.659458e-6;
    let est_time: Time = Time::new(4.0*time_per_spp*(samples*width*height) as f64);

    println!("Estimated time [ DEBUG ]: {}", est_time.get_time());

    let time_per_spp: f64 = 0.432191e-6;
    let est_time: Time = Time::new(4.0*time_per_spp*(samples*width*height) as f64);

    println!("Estimated time [RELEASE]: {}", est_time.get_time());


    let threadpool = ThreadPool::new(num_threads);
    let (tx, rx) = channel();

    let scene = Arc::new(Scene::new2(width, height));

    let time_start = time::precise_time_s();
    for y in 0..height {
        let tx = tx.clone();
        let scene = scene.clone();

        threadpool.execute(move || {
            let line = single_row(y, width, height, samples, scene);
            tx.send((y, line)).unwrap();
        });
    }

    let mut left = height;
    let mut screen : Vec<Vec<Vec3d>> = Vec::new();
    for _y in 0..height {
        screen.push(Vec::new());
    }
    while left > 0 {
        print!("Rendering ({} spp) {:.4}%...\r", samples * 4, 100.0 * (height - left) as f64 / height as f64);
        io::stdout().flush().ok().expect("Could not flush stdout");
        let (y, line) = rx.recv().unwrap();
        screen[y] = line;
        left -= 1;
    }

    let time_taken = time::precise_time_s() - time_start;
    println!("Finished rendering. Time taken: {}", Time::new(time_taken).get_time());
    println!("DEBUG time_per_spp: {}", (time_taken as f64/(width*height*4*samples) as f64)*1e6);

    let image_name = format!("image_{}_{}_{}.png", width, height, samples*4);
    save_image(&screen, &image_name);
}


pub fn save_image(image: &Vec<Vec<Vec3d>>, image_name: &str) {
    let height = image.len();
    let width = image[0].len();

    let mut buffer = Vec::<u8>::with_capacity((width*height*3) as usize);
    unsafe { buffer.set_len((width*height*3) as usize); }


    for h in 0..height {
        for w in 0..width {
            let h2 = height - h - 1;
            buffer[h2*width*3 + w*3] = to_u8(image[h][w].x);
            buffer[h2*width*3 + w*3 + 1] = to_u8(image[h][w].y);
            buffer[h2*width*3 + w*3 + 2] = to_u8(image[h][w].z);
        }
    }

    image::save_buffer(&Path::new(image_name), &buffer, width as u32, height as u32, image::RGB(8))
        .ok().expect("Failed to save the image");
}




#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_row() {
        let width = 100;
        let height = 100;
        let samples = 10;
        let scene = Arc::new(Scene::new2(width, height));
        let result = single_row(0, width, height, samples, scene.clone());
        assert_eq!(result.len(), width);
    }

/*
    #[test]
    fn test_intersection() {
        let sphere = Sphere {radius:1e5, position: Vec3d{x:1e5+1.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.25,z:0.25}, reflection: ReflectType::DIFF};
        let ray = Ray { origin: Vec3d { x: 14.234, y: 46.039, z: 14.234 }, direction: Vec3d { x: -0.702, y: -0.117, z: -0.702 } };

        //assert_eq!(sphere.intersect(&ray), 18.85);
    }

    #[test]
    fn test_light() {
        let scene = Scene::new();
        let ray = Ray { origin: Vec3d { x: 50.0, y: 50.0, z: 81.6 }, direction: Vec3d { x: -0.0, y: 1.0, z: 0.0 } }; // aim directly at light
        //spheres.push(Sphere {radius:600.0, position: Vec3d{x:50.0,y:681.6-0.27,z:81.6}, emission: Vec3d{x:12.0,y:12.0,z:12.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light

        let (intersected, dist, id) = scene.intersect(&ray);
        assert!(float_eq(dist, 31.33));

        let result = scene.radiance(&ray, 0, 1);
        assert!(result != Vec3d::zeros());
    }*/
}

#[cfg(all(feature = "unstable", test))]
mod bench {
    extern crate test;

    use super::*;
    use self::test::Bencher;

    #[bench]
    fn bench_single_row(b: &mut Bencher) {
        let width = 100;
        let height = 100;
        let samples = 10;
        let scene = Arc::new(Scene::new2(width, height));
        b.iter(|| single_row(5, width, height, samples, scene.clone()));
    }

    #[bench]
    fn bench_update_pixel(b: &mut Bencher) {
        let width = 100;
        let height = 100;
        let samples = 10;
        let scene = Arc::new(Scene::new2(width, height));
        b.iter(|| update_pixel(5.0, 5.0, width as f64, height as f64, samples, &scene));
    }

}
