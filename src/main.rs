extern crate rand;
extern crate threadpool;
extern crate num_cpus;
extern crate clock_ticks;
extern crate docopt;

use std::env;
use std::f64;
use std::io::BufWriter;
use std::fs::{OpenOptions};
use std::io::Write;
use std::io;
use std::sync::Arc;
use threadpool::ThreadPool;
use std::sync::mpsc::channel;
use docopt::Docopt;

mod vector_3d;

use vector_3d::Vec3d;

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
	} else 	if x > 1.0 {
		return 1.0;
	}
	x
}

pub fn to_int(x: f64) -> i32 {
	(clamp(x).powf(1.0/2.2)*255.0+0.5) as i32
}


pub struct Scene {
	pub spheres: Vec<Sphere>
}

impl Scene {
	pub fn new() -> Scene {
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
			spheres: spheres
		}
	}

	pub fn new2() -> Scene {
		let mut spheres = Vec::new();
		spheres.push(Sphere {radius:1e5, position: Vec3d{x:1e5+1.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.8,y:0.1,z:0.1}, reflection: ReflectType::DIFF}); // left
		spheres.push(Sphere {radius:1e5, position: Vec3d{x:-1e5+99.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.1,y:0.3,z:0.70}, reflection: ReflectType::DIFF}); // right
		spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:40.8,z:1e5}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // back
		spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:40.8,z:-1e5+170.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // front
		spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:1e5,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // bottom
		spheres.push(Sphere {radius:1e5, position: Vec3d{x:50.0,y:-1e5+81.6,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.75,z:0.75}, reflection: ReflectType::DIFF}); // Top
		spheres.push(Sphere {radius:16.5, position: Vec3d{x:27.0,y:16.5,z:47.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.95, reflection: ReflectType::SPEC}); // mirror
		spheres.push(Sphere {radius:16.5, position: Vec3d{x:73.0,y:30.0,z:78.0}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:1.0,y:1.0,z:1.0}*0.95, reflection: ReflectType::REFR}); // glass
		spheres.push(Sphere {radius:600.0, position: Vec3d{x:50.0,y:681.6-0.27,z:81.6}, emission: Vec3d{x:12.0,y:12.0,z:12.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light

		Scene {
			spheres: spheres
		}
	}

	pub fn intersect(&self, ray: &Ray) -> (bool, f64, i32) {
		let n = self.spheres.len();
		let mut d: f64;
		let mut t: f64 = 1e20;
		let mut id: i32 = 0;
		let inf = 1e20;

		for i in (0..n).rev() {
			//println!("i: {}",i);
			d = self.spheres[i].intersect(ray);
			if d != 0.0 && d < t {
				//println!("ray intersected sphere {}: {:?}", i, ray);
				t = d;
				id = i as i32;
			}
		}
		//println!("intersects - t: {}, inf: {}", t, inf);
		(t < inf, t, id)
	}

	pub fn radiance(&self, ray: &Ray, mut depth: i32) -> Vec3d {
		//println!("radiance - depth: {},  ray: {:?}", depth, ray);
		let (intersects, t, id) = self.intersect(ray);
		if !intersects {  // if miss, return black
			//println!("No Intersect!");
			return Vec3d{x:0.0, y:0.0, z:0.0};
		}

		let obj: Sphere = self.spheres[id as usize].clone();  // the hit object

		let x: Vec3d = ray.origin + ray.direction*t; // point on sphere where intersects
		let n: Vec3d = (&x - &obj.position).normalise(); // surface normal of intersection point
		let nl: Vec3d; // corrected normal (ie internal or external intersection)
		if n.dot(&ray.direction) < 0.0 { // dot product negative if ray is internal
			nl = n;
		} else {
			nl = n*-1.0;
		}

		//println!("ray: {:?}",ray);
		//println!("t: {:?}",t);
		//println!("obj: {:?}",obj);
		//println!("x: {:?}",x);
		//println!("n: {:?}",n);
		//println!("nl: {:?}",nl);

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
		if depth > 5 {
			if rand::random::<f64>()+0.5 < p {
				f = f*(1.0/p); // normalise colour [0,1]
			} else {
				// This might be at the wrong if statement
				return obj.emission;
			}
		}

		match obj.reflection {
			ReflectType::DIFF => {
				let r1: f64 = 2.0*f64::consts::PI*rand::random::<f64>();
				let r2: f64 = rand::random();
				let r2s = r2.sqrt();

				let w: Vec3d = nl.clone();
				let u: Vec3d;
				if w.x.abs() > 0.1 {
					u = (Vec3d{x:0.0, y:1.0,z:0.0}.cross(w)).normalise();
				} else {
					u = (Vec3d{x:1.0, y:0.0,z:0.0}.cross(w)).normalise();
				}
				//println!("obj.emission DIFF w: {:?}, u: {:?}", w, u);
				assert!((u.length() - 1.0).abs() < 0.001);

				let v: Vec3d = w.cross(u);

				let d: Vec3d = (u*r1.cos()*r2s + v*r1.sin()*r2s + w*(1.0-r2).sqrt()).normalise();
				assert!((d.length() - 1.0).abs() < 0.001);
				//println!("obj.emission DIFF u: {:?}, v: {:?}, w: {:?}, r1: {:?}, r2s: {:?}", u, v, w, r1, r2s);
				let fmu = f.mult(self.radiance(&Ray{origin: x, direction: d}, depth));
				//println!("obj.emission DIFF - fmu.x: {}, fmu.y: {}, fmu.z: {}, depth: {}", fmu.x, fmu.y, fmu.z, depth);
				//println!("obj.emission.x: {}, obj.emission.y: {}, obj.emission.z: {}", obj.emission.x, obj.emission.y, obj.emission.z);
				return obj.emission + fmu;
			},
			ReflectType::SPEC => {
				//println!("obj.emission SPEC, depth: {}", depth);
				let d: Vec3d = ray.direction-(n*2.0*n.dot(&ray.direction));
				return obj.emission + f.mult(self.radiance(&Ray{origin: x, direction: d}, depth));
			},
				_ => {}

		}

		let refl_ray: Ray = Ray{origin: x, direction: ray.direction-n*2.0*n.dot(&ray.direction)}; // Ideal dielectric REFRACTION
		let into: bool = n.dot(&nl) > 0.0;

		let nc: f64 = 1.0;
		let nt: f64 = 1.5;
		let nnt: f64;
		if into { // Ray from outside going in?
			nnt = nc/nt;
		} else {
			nnt = nt/nc;
		}
		let ddn: f64 =ray.direction.dot(&nl);
		let cos2t: f64 = 1.0 - nnt * nnt * (1.0 - ddn * ddn);
		if cos2t < 0.0 { // Total internal reflection
			//println!("obj.emission cos2t < 0.0");
			return obj.emission + f.mult(self.radiance(&refl_ray, depth));
		}

		let tdir: Vec3d;

		if into {
			tdir = (ray.direction*nnt - n*(1.0*(ddn*nnt+cos2t.sqrt()))).normalise();
		} else {
			tdir = (ray.direction*nnt - n*(-1.0*(ddn*nnt+cos2t.sqrt()))).normalise();
		}

		let a: f64 = nt- nc;
		let b: f64 = nt + nc;
		let r0: f64 = a*a/(b*b);
		let c: f64;

		if into {
			c = 1.0 + ddn;
		} else {
			c = 1.0 - tdir.dot(&n);
		}

		let re: f64 = r0+(1.0-r0)*c.powi(5);
		let tr: f64 = 1.0 - re;
		let p: f64 = 0.25+0.5*re;
		let rp: f64 = re/p;
		let tp: f64 = tr/(1.0-p);

		if depth > 2 {
			if rand::random::<f64>() < p {
				//println!("rand::random::<f64>() < p - depth: {}", depth);
				return obj.emission + f.mult(self.radiance(&refl_ray, depth)*rp);
			} else {
				//println!("rand::random::<f64>() >= p - depth: {}", depth);
				return obj.emission + f.mult(self.radiance(&Ray{origin: x, direction: tdir}, depth)*tp);
			}
		} else {
			//println!("reflect ray: {}", depth);
			return obj.emission + f.mult(self.radiance(&refl_ray, depth)*re+self.radiance(&Ray{origin: x, direction: tdir}, depth)*tr);
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
	let samps: usize = args.get_str("<samples>").parse().unwrap_or(100);

	println!("width: {}, height: {}, samples: {}", width, height, samps);


	let num_threads = num_cpus::get();

	let time_per_spp: f64 = 3.659458e-6;
	let est_time: Time = Time::new(4.0*time_per_spp*(samps*width*height) as f64);

	println!("Estimated time [ DEBUG ]: {}", est_time.get_time());

	let time_per_spp: f64 = 0.432191e-6;
	let est_time: Time = Time::new(4.0*time_per_spp*(samps*width*height) as f64);

	println!("Estimated time [RELEASE]: {}", est_time.get_time());

	let scene = Arc::new(Scene::new2());
	let cam: Ray = Ray{origin: Vec3d{x:50.0,y:50.0,z:295.6}, direction: Vec3d{x:0.0,y:-0.042612, z:-1.0}.normalise()};
	//let cam: Ray = Ray{origin: Vec3d{x:50.0,y:42.0,z:235.6}, direction: Vec3d{x:0.0,y:-0.042612, z:-1.0}.normalise()};

	let fov: f64 = 0.5135;
	let cx: Vec3d = Vec3d{x:(width as f64)*fov/(height as f64),y:0.0,z:0.0}; // x direction increment
	let cy: Vec3d = (cx % cam.direction).normalise()*fov;                    // y direction increment

	let threadpool = ThreadPool::new(num_threads);
	let (tx, rx) = channel();

	let mut image = Vec::<Vec3d>::with_capacity((width*height) as usize);
	unsafe { image.set_len((width*height) as usize); }

	let time_start = clock_ticks::precise_time_s();

	for y in 0..height {
		let tx = tx.clone();
		let scene = scene.clone();
		let mut r: Vec3d = Vec3d{x:0.0, y: 0.0, z: 0.0};

		threadpool.execute(move || {
			let mut line = Vec::with_capacity(width);
			for x in 0..width {
				let mut sum = Vec3d{x:0.0,y:0.0,z:0.0};
				for sy in 0..2 {
					for sx in 0..2 {
						for _ in 0..samps {
							let r1: f64 = 2.0*rand::random::<f64>(); //erand48(xi);
							let r2: f64 = 2.0*rand::random::<f64>(); //erand48(xi);

							let dx: f64;
							let dy: f64;

							if r1 < 1.0 {
								dx = r1.sqrt() - 1.0;
							} else {
								dx = 1.0 - (2.0-r1).sqrt();
							}

							if r2 < 1.0 {
								dy = r2.sqrt() - 1.0;
							} else {
								dy = 1.0 - (2.0-r2).sqrt();
							}

							let mut d: Vec3d = cx*((((sx as f64)+0.5 + dx)/2.0 + (x as f64)) / (width as f64) - 0.5) +
											cy*((((sy as f64)+0.5 + dy)/2.0 + (y as f64)) / (height as f64) - 0.5) + cam.direction;
							//println!("original dir: {:?}", d);
							let rad: Vec3d = scene.radiance(&Ray{origin: cam.origin + d*140.0, direction: d.normalise()},0);
							//println!("rad.x: {}, rad.y: {}, rad.z: {}", rad.x, rad.y, rad.z);
							r = r + rad*(1.0/samps as f64);

						}
						let v: Vec3d = Vec3d{x: clamp(r.x), y: clamp(r.y), z: clamp(r.z)};
						//println!("v.x: {}, v.y: {}, v.z: {}", v.x, v.y, v.z);
						sum = sum + v*0.25;
						r = Vec3d{x:0.0,y:0.0,z:0.0};
					}
				}
				line.push(sum);
			}
			tx.send((y, line)).unwrap();
		});
	}



	let mut left = height;
	let mut screen : Vec<Vec<Vec3d>> = Vec::new();
	for _y in 0..height {
		screen.push(Vec::new());
	}
	while left > 0 {
		print!("Rendering ({} spp) {:.4}%...\r", samps * 4, 100.0 * (height - left) as f64 / height as f64);
		io::stdout().flush().ok().expect("Could not flush stdout");
		let (y, line) = rx.recv().unwrap();
		screen[y] = line;
		left -= 1;
	}

	let time_taken = clock_ticks::precise_time_s() - time_start;
	println!("Finished rendering. Time taken: {}", Time::new(time_taken).get_time());
	println!("DEBUG time_per_spp: {}", (time_taken as f64/(width*height*4*samps) as f64)*1e6);

	// We create file options to write
	let file = OpenOptions::new().write(true).create(true).open(format!("image_{}_{}_{}.ppm", width, height, samps*4)).unwrap();

	// We create a buffered writer from the file we get
	let mut writer = BufWriter::new(&file);
	// Then we write to the file. write_all() calls flush() after the write as well.
	let mut b = format!("P3\n{} {}\n{}\n", width, height, 255).into_bytes();
	writer.write_all(&b);
	for i in (0..height as usize).rev() {
		for j in 0..width as usize {
			b = format!("{} {} {}\n", to_int(screen[i][j].x), to_int(screen[i][j].y), to_int(screen[i][j].z)).into_bytes();
			writer.write_all(&b);
		}
	}
}





#[cfg(test)]
mod test {
	use super::*;
	use vector_3d::Vec3d;
	use rand;

	#[test]
	fn test_intersection() {
		let sphere = Sphere {radius:1e5, position: Vec3d{x:1e5+1.0,y:40.8,z:81.6}, emission: Vec3d{x:0.0,y:0.0,z:0.0}, color: Vec3d{x:0.75,y:0.25,z:0.25}, reflection: ReflectType::DIFF};
		let ray = Ray { origin: Vec3d { x: 14.234, y: 46.039, z: 14.234 }, direction: Vec3d { x: -0.702, y: -0.117, z: -0.702 } };

		println!("sphere: {:?}", sphere);
		println!("ray: {:?}", ray);
		println!("sphere.intersect(&ray): {:?}", sphere.intersect(&ray));
		//assert_eq!(sphere.intersect(&ray), 18.85);
	}

	#[test]
	fn test_light() {
		let scene = Scene::new();
		let ray = Ray { origin: Vec3d { x: 50.0, y: 50.0, z: 81.6 }, direction: Vec3d { x: -0.0, y: 1.0, z: 0.0 } }; // aim directly at light
		//spheres.push(Sphere {radius:600.0, position: Vec3d{x:50.0,y:681.6-0.27,z:81.6}, emission: Vec3d{x:12.0,y:12.0,z:12.0}, color: Vec3d{x:0.0,y:0.0,z:0.0}, reflection: ReflectType::DIFF}); // light

		let (intersected, dist, id) = scene.intersect(&ray);
		assert!(float_eq(dist, 31.33));

		let result = scene.radiance(&ray, 0);
		println!("result: {:?}", result);
		assert!(result != Vec3d::zeros());
	}

	#[test]
	fn test_light_scene() {
		let x = 30;
		let y = 43;

		let width = 50;
		let height = 50;
		let samps = 1;

		let scene = Scene::new();
		let mut r: Vec3d = Vec3d{x:0.0, y: 0.0, z: 0.0};

		let cam: Ray = Ray{origin: Vec3d{x:50.0,y:50.0,z:295.6}, direction: Vec3d{x:0.0,y:-0.042612, z:-1.0}.normalise()};

		let cx: Vec3d = Vec3d{x:(width as f64)*0.5135/(height as f64),y:0.0,z:0.0}; // x direction increment
		let cy: Vec3d = (cx % cam.direction).normalise()*0.5135;                    // y direction increment

		let mut sum = Vec3d{x:0.0,y:0.0,z:0.0};
		for sy in 0..2 {
			for sx in 0..2 {
				for _ in 0..samps {
					let r1: f64 = 2.0*rand::random::<f64>(); //erand48(xi);
					let r2: f64 = 2.0*rand::random::<f64>(); //erand48(xi);

					let dx: f64;
					let dy: f64;

					if r1 < 1.0 {
						dx = r1.sqrt() - 1.0;
					} else {
						dx = 1.0 - (2.0-r1).sqrt();
					}

					if r2 < 1.0 {
						dy = r2.sqrt() - 1.0;
					} else {
						dy = 1.0 - (2.0-r2).sqrt();
					}

					let mut d: Vec3d = cx*((((sx as f64)+0.5 + dx)/2.0 + (x as f64)) / (width as f64) - 0.5) +
					cy*((((sy as f64)+0.5 + dy)/2.0 + (y as f64)) / (height as f64) - 0.5) + cam.direction;
					println!("original dir: {:?}", d);
					let rad: Vec3d = scene.radiance(&Ray{origin: cam.origin + d*140.0, direction: d.normalise()},0);
					println!("rad: {:?}", rad);
					r = r + rad*(1.0/samps as f64);
					println!("r: {:?}", r);

				}
				let v: Vec3d = Vec3d{x: clamp(r.x), y: clamp(r.y), z: clamp(r.z)};
				println!("v: {:?}", v);
				sum = sum + v*0.25;
				println!("sum: {:?}", sum);
				r = Vec3d{x:0.0,y:0.0,z:0.0};
			}
		}

		assert!(sum != Vec3d::zeros());
	}
}
