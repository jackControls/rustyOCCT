use rusty_occt::{Boundary, Frame3, Point2, Profile, Solid, Tolerance};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tolerance = Tolerance::default();
    let outer = Boundary::rectangle(80.0, 40.0, tolerance)?;
    let holes = [(10.0, 10.0), (70.0, 10.0), (70.0, 30.0), (10.0, 30.0)]
        .into_iter()
        .map(|(x, y)| Boundary::circle(Point2::new(x, y), 3.0, tolerance))
        .collect::<Result<Vec<_>, _>>()?;
    let profile = Profile::new(outer, holes, tolerance)?;
    let plate = Solid::extrude(profile, Frame3::xy(), 0.0, 6.0)?;
    let mass = plate.mass_properties();
    println!("Plate: 80 x 40 x 6 mm; four exact 6 mm through holes");
    println!("Volume: {:.6} mm^3", mass.volume);
    println!("Surface area: {:.6} mm^2", mass.surface_area);
    println!("Centroid: {:?} mm", mass.centroid.to_array());
    println!(
        "Topology: {} vertices, {} edges, {} faces",
        plate.topology().vertices().len(),
        plate.topology().edges().len(),
        plate.topology().faces().len()
    );
    println!(
        "Closed shell verified; Euler characteristic {}",
        plate.topology().euler_characteristic()
    );
    Ok(())
}
