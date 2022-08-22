"""
Makes slab-ligand combinations for Sn based stuff.
    Usage:
        Initialize and object.
        Set the required attributes for the ligands
        use prepare methods to prepare the acenes after setting required settings
        use the make ligands to make the ligand
"""

import ase

class Slab_Sn_Ligand():
    def __init__(self) -> None:
        self.XC = "LDA"
        self.number_of_repetitions_in_ligand = 3

        # Setting the lattice constant even though the exact value is not needed
        if self.XC == "LDA": a0 = 6.479  # This is the lattice constant for the slab

    def prepare_acenes(self) -> None:
        pass

    def prepare_PPV(self) -> None:
        pass

    def make_structure(self) -> ase.Atoms:
        a0 = 6.479  # This is the lattice constant for the slab
        print(f"now working on ligand:{ligand}")

        relaxed_structure_UC = UC-1
        if passivation: passivation_text = f"H"
        else: passivation_text = ""
        if direction != "111": terminating_surface = ""

        if method == "R":
            if relaxed_name == "":
                # Trying to automatically get the relaxed name.
                relaxed_name = f"Sn{direction}{terminating_surface}_{relaxed_structure_UC}{passivation_text}"
            print("UC and relaxed_structure_UC:", UC, relaxed_structure_UC)
            relaxed_structure_atoms = read(f"{base_folder}/{relaxed_name}^Calc=SIESTA^Struct=bulk^Specie=Sn-^KE=20^K=2^R=80^a=6.47^type=scf/{relaxed_name}.xyz")
            relaxed_slab = MakeSlab(a0 = a0)
            relaxed_slab.load(relaxed_structure_atoms)
            # relaxed_structure_atoms.edit()

        # Setting the label
        # THis might be unecessary since the terminating surface is already set to null if surface is not 111
        # Wait for working code for unrelaxed runs and then adjust (remove and shorten the code here)
        if direction == "111":
            label = f"Sn{direction}{terminating_surface}_{int(UC)}"
        else:
            label = f"Sn{direction}_{int(UC)}"
        
        # Label but we also set the species here and depending on how you made it, it might be different.
        if passivation:
            label = f"{label}H"
        if attach_ligand:
            label = f"{label}{ligand}"

        if method != "R":
            # This was when we used to make the full thickness together
            if direction == "100":
                slab = Diamond100(6.479, UC*1)
            if direction == "110":
                slab = Diamond110(6.479, UC*1)
            if direction == "111":
                slab = Diamond111(6.479, UC*1)
                # slab = Diamond111_diamond(6.479, UC*1)
        else:
            number_of_UCs_required = UC
            if direction == "100":
                slab = Diamond100(6.479, number_of_UCs_required)
            if direction == "110":
                slab = Diamond110(6.479, number_of_UCs_required)
            if direction == "111":
                slab = Diamond111(6.479, number_of_UCs_required)
                # slab = Diamond111_diamond(6.479, UC*1)
            if direction == "210":
                slab = Diamond210(6.479, 1)

        # if method == "R":
            # this method is not used now since it was really hard to do
            # we are now useing the update_to_relaxed technique that you will find further down in the code
        #     relaxed_slab.sandwich_slab_with_bread(slab.atoms, bread_y_offset=a0)
        #     # relaxed_slab.sandwich_slab_with_bread(slab.atoms)
        #     write(f"split_structure.vasp", relaxed_slab.atoms)
        #     find_cell(relaxed_slab.atoms, a0)
        #     relaxed_slab.center_to_cell()
        #     # relaxed_slab.edit()
        #     relaxed_slab.make_inversion_symmetric()
        #     # relaxed_slab.edit()
        #     # LI = relaxed_slab
        # slab.atoms.edit()

        if method != "R":
            # Here we will make is so that we can decide if the surface of the 111 direction is terminated as A or B
            if direction == "111" and terminating_surface == "B":
                largest_index = len(slab.atoms) - 1
                del slab.atoms[largest_index - 4]
                del slab.atoms[largest_index - 6]
                if UC == 1: del slab.atoms[6]
                else: del slab.atoms[8]
                del slab.atoms[0]
                # slab.atoms.edit()
            
        # slab.atoms.edit()
        number_of_Sn_atoms = len(slab.atoms)
        # print(number_of_Sn_atoms)
        atoms_to_constrain = list(range(1,number_of_Sn_atoms + 1))

        # if passivation and method != "R":
        if passivation:
            if direction == "111" and terminating_surface == "B":
                slab.passivate_zinc_blende_slab(1.7, ["+z"], slab_miller_index=direction, slab_termination="B")
            else:
                slab.passivate_zinc_blende_slab(1.7, ["+z"], slab_miller_index=direction)
        
        if not attach_ligand and not method == "R":
            # We have to do this at the end so that we can make sure that we are flipping the final structure
            # This is for the part when you do not need to attach a ligand
            if cell_inversion_symmetric:
                # slab.edit()
                # slab.center_to_cell()
                slab.add_vacuum(5)
                # slab.edit()
                slab.make_inversion_symmetric()
                # slab.edit()
                kwargs = {"":0}
                if direction == "111":
                    kwargs.update({"x_offset":-a0/4})
                    kwargs.update({"y_offset":-a0/4})
                    kwargs.update({"x":-a0/2})
                    kwargs.update({"x_offset":0})
                    kwargs.update({"y_offset":0})
                if direction == "100":
                    kwargs.update({"x_offset":a0/4})
                    kwargs.update({"y_offset":-a0/4})
                slab.set_cell(**kwargs)
                slab.center_to_cell()
            else:
                # slab.add_vacuum(11)
                slab.center_to_cell()
            # slab.atoms.edit()
            slab.add_vacuum(11)
            LI = slab.copy()


        if attach_ligand:
            # Here we have the part where we attach ligands
            # print(ligand)
            LI = Insert_ligand(name=ligand)
            # LI.orient([0,0,-1])  # Orienting the ligand along the z axis
            # write(f"{ligand}_aligned_to_z.xyz", LI.atoms)
            # LI.edit()
            # LI.atoms.rotate(-90,[0,0,1],(0,0,0))
            # LI.atoms.rotate(10,[1,0,0],(0,0,0))
            LI.atoms.rotate(-90,[0,1,0],(0,0,0))
            # LI.edit()
            if "Ac" in ligand:
                LI.atoms.rotate(20,[1,0,0],(0,0,0))
                LI.atoms.rotate(90,[0,0,1],(0,0,0))
                LI.atoms.rotate(10,[0,1,0],(0,0,0))
                LI.atoms.rotate(-85,[0,0,1],(0,0,0))
                # LI.atoms.rotate(15,[0,1,0],(0,0,0))
                # LI.atoms.rotate(10,[0,0,1],(0,0,0))
            # LI.edit()
            # LI.atoms.rotate(25,[0,1,0],(0,0,0))
            # LI.edit()
            if not "Ac" in ligand:
                # LI.edit()
                LI.align_atoms_along_axis(sites[ligand][1], sites[ligand][2], [0,0,0], [0,0,1])  #  Here we are aligning the two S atoms
                # LI.edit()
            write(f"{ligand}_rotated_to_z_axis.xyz", LI.atoms)

            # slab.atoms.edit()
            # attach_through = ligand_attach_through[ligand]
            # atoms_to_delete = ligand_atoms_to_delete[ligand]
            if direction == "100":
                slab_site_bulk = len(slab)-5
            elif direction == "110":
                slab_site_bulk = len(slab)-5
            # print(f"Feeding in sb: {slab_site_bulk}")
            slab_site_passivant = slab_site_bulk + 3
            LI.attach_to_slab_Qunfei_dense(slab.atoms, slab_site_bulk, slab_site_passivant, sites[ligand][0], sites[ligand][3], sites[ligand][1], sites[ligand][2])
            # LI.edit()
            # LI.align_atoms_along_axis(sites[ligand][1], sites[ligand][2], [0,0,0], [0,0,1])  #  Here we are aligning the two S atoms
            # LI.edit()
            LI.make_inversion_symmetric()
            # LI.edit()
            kwargs = {}
            kwargs.update({"x_offset":a0/4})
            kwargs.update({"y_offset":a0/4})
            find_cell(LI.atoms, a0, **kwargs)
            if method == "R":
                LI.update_to_relaxed_coordinates(relaxed_slab.atoms)
            # LI.edit()
            extreme_coordinates = get_extreme_coordinates(LI.atoms)
            LI = LI.atoms # this is because if run is of only slabs and no ligands the object is a atoms type object
            # print("cell:" , LI.get_cell())