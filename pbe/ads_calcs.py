#!/usr/bin/env python3

import click
import yaml


@click.group()
def cli():
    """Tools to process data for metal-substrate calculations."""
    pass


@click.command()
@click.option(
    "-f",
    "--filename",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    help="path to geometry file",
)
@click.option(
    "-b",
    "--bulk_atom",
    type=str,
    required=True,
    help="element of the atoms in the bulk",
)
@click.option(
    "-c",
    "--constr_layers",
    "n_constr_layers",
    type=int,
    required=True,
    help="number of constrained layers in the bulk",
)
def reorder_atoms(filename, bulk_atom, n_constr_layers):
    """Sort bulk atoms by ascending z-coordinate."""

    with open(filename, "r") as geom:
        lines = geom.readlines()

    init_at_line = 0

    for i in range(len(lines)):
        if "atom" == lines[i].split()[0]:
            init_at_line = i
            break

    z_lines = []

    # Find all the z coordinates and sort them from smallest to largest
    for line_content in lines:
        if bulk_atom in line_content.split():
            z_lines.append((float(line_content.split()[3]), line_content))

    z_lines.sort(key=lambda elem: elem[0])

    # Get non-repeated list of z coors of bulk atoms
    largest_coor = 0
    coors = [0.0000000000000000]
    for tup in z_lines:
        coor = tup[0]

        if largest_coor < coor:
            coors.append(coor)
            largest_coor = coor

    # Save the bottom n layers
    constr_layers = coors[:n_constr_layers]

    # Insert the new atoms and re-apply the constraints
    insert_counter = 0
    for i, tup in zip(range(init_at_line, len(z_lines) + init_at_line), z_lines):
        lines[i + insert_counter] = tup[1]

        if float(lines[i + insert_counter].split()[3]) in constr_layers:
            insert_counter += 1
            lines[i + insert_counter] = "    constrain_relaxation .true.\n"

    with open(filename, "w") as geom:
        geom.writelines(lines)


@click.command()
@click.option(
    "-i",
    "--inputfile",
    type=click.Path(exists=True, file_okay=True, dir_okay=False),
    required=True,
    help="path to geometry file",
)
@click.option(
    "-b",
    "--bulk_atom",
    type=str,
    required=True,
    help="element of the atoms in the bulk",
)
@click.option(
    "-p",
    "--ads_per_atom",
    default=False,
    show_default=True,
    is_flag=True,
    help="print the adsorption height per atom",
)
@click.option(
    "-a",
    "--substrate_atoms",
    "ads_atoms",
    multiple=True,
    type=str,
    required=True,
    help="specify each atom in the adsorbate",
)
@click.option(
    "-l",
    "--bulk_layers",
    type=int,
    required=True,
    help="number of layers in the bulk",
)
def substrate_info(inputfile, bulk_atom, bulk_layers, ads_per_atom, ads_atoms):
    """Calculate substrate information."""

    with open(inputfile, "r") as geom:
        lines = geom.readlines()

    ads_z_coors = []
    bulk_z_coors = []
    tot_ad_atoms_ref = []

    # Find the z-coordinates of each atom
    for atom in ads_atoms:
        n_atoms = []
        ad_atoms_ref = []
        for line in lines:
            spl = line.split()
            if atom == spl[-1]:
                n_atoms.append(float(spl[3]))
                ad_atoms_ref.append(spl[-1])

            if bulk_atom == spl[-1] and spl[-1] not in bulk_z_coors:
                bulk_z_coors.append(float(line.split()[3]))

        ads_z_coors.append(n_atoms)
        tot_ad_atoms_ref.append(ad_atoms_ref)

    # Find the z-coordinate of the 2nd layer of bulk atoms
    nd_bulk_layers = list(dict.fromkeys(bulk_z_coors))
    hyp_bulk_height = nd_bulk_layers[1] * (bulk_layers - 1)

    if ads_per_atom is True:
        print()
        print(" Atom | Adsorption Height (Å) ")
        print(" -----|----------------------")

    ads_heights = []
    heighest_atom = 0
    lowest_atom = max([i for i in ads_z_coors])

    # Iterate over the z-coordinates of each adsorbate atom
    for i, element in enumerate(ads_z_coors):
        for j, coor in enumerate(element):
            ads_height = ads_z_coors[i][j] - hyp_bulk_height

            # Calculate the hyp bulk height and adsorption height per atom
            if ads_per_atom is True:
                print(" ", tot_ad_atoms_ref[i][j], "  |", ads_height)

                if element == ads_z_coors[-1] and coor == element[-1]:
                    print()

            ads_heights.append(ads_height)

    print("Hypothetical unrelaxed bulk height:", round(hyp_bulk_height, 4), "Å")

    # Find the A_r of each specified adsorbate atom
    with open("/home/dylanmorgan/comp_chem/warwick/azupyrene_cu/elements.yml", "r") as elements:
        elements = yaml.load(elements, Loader=yaml.FullLoader)

    # Calculate the weighteed average adsorption heights by A_r
    weighted_heights = []
    at_masses = []
    for height, el in zip(ads_heights, ads_atoms):
        A_r = float(elements[el]["atomic_mass"])
        at_masses.append(A_r)
        weighted_heights.append(A_r * height)

    # Print the weighted average adsorption height of all adsorbate atoms
    wavg_ads_height = sum(weighted_heights) / sum(at_masses)
    print("Weighted average atomic adsorption height:", round(wavg_ads_height, 4), "Å")

    # Print the weighted average adsorpation height from the hypothetical bulk height
    ads_height_nixsw = hyp_bulk_height - wavg_ads_height
    print(
        "Theoretically measured weighted average substrate height:",
        round(ads_height_nixsw, 4),
        "Å",
    )

    # Take the difference of the heighest and lowest adsorbate atom
    heighest_atom = max([max(sublist) for sublist in ads_z_coors])
    lowest_atom = min([min(sublist) for sublist in ads_z_coors])

    delta_ads = heighest_atom - lowest_atom
    print("Substrate strain:", round(delta_ads, 4), "Å")


@click.command()
@click.option(
    "-f",
    "--full_geom_opt",
    "full_go",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="path to directory of full geometry optimisation",
)
@click.option(
    "-so",
    "--surface_opt",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="path to directory of optimised surface reference optimisation",
)
@click.option(
    "-ss",
    "--surface_sp",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="path to directory of single-point surface reference optimisation",
)
@click.option(
    "-ao",
    "--substrate_opt",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="path to directory of optimised substrate reference optimisation",
)
@click.option(
    "-as",
    "--substrate_sp",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="path to directory of single-point substrate reference optimisation",
)
@click.option(
    "-a",
    "--substrate_atoms",
    "ads_atoms",
    multiple=True,
    type=str,
    required=True,
    help="specify each atom in the adsorbate",
)
def parse_energies(
    full_go, surface_opt, surface_sp, substrate_opt, substrate_sp, ads_atoms
):
    """Parse energies from the various geometry optimisations."""

    parse_line = "| Final zero-broadening corrected energy (caution - metals only)"
    outfiles = (
        f"{full_go}/aims.out",
        f"{surface_opt}/aims.out",
        f"{substrate_opt}/aims.out",
        f"{surface_sp}/aims.out",
        f"{substrate_sp}/aims.out",
    )
    energies = [0.0, 0.0, 0.0, 0.0, 0.0]

    for fn, file in enumerate(outfiles):
        with open(file, "r") as out:
            out_lines = out.readlines()

        for line in out_lines:
            spl = line.split()
            if parse_line in line:
                energies[fn] = float(spl[-2])

    # Find the total number of adsorbate atoms
    with open(f"{substrate_opt}/geometry.in", "r") as geom:
        lines = geom.readlines()

    n_ads_atoms = []
    for ads in ads_atoms:
        ads_atom_count = 0
        for line in lines:
            spl = line.split()
            if "atom" in spl and ads in spl:
                ads_atom_count += 1

        n_ads_atoms.append(ads_atom_count)

    # Calculate the adsorption energy
    ads_e = energies[0] - energies[1] - energies[2]
    ads_i = energies[0] - energies[3] - energies[4]

    print("Full geometry optimisation energy:", round(energies[0], 4), "eV")
    print("Surface energy:", round(energies[1], 4), "eV")
    print("Substrate energy:", round(energies[2], 4), "eV")
    print("Adsorption energy:", round(ads_e, 4), "eV")

    for n, ads in zip(n_ads_atoms, ads_atoms):
        print(
            f"Adsorption energy per substrate {ads}:",
            round(1000 * (ads_e / n), 4),
            "meV",
        )

    print("Interaction energy:", round(ads_i, 4), "eV")


if __name__ == "__main__":
    cli.add_command(reorder_atoms)
    cli.add_command(substrate_info)
    cli.add_command(parse_energies)
    cli()
