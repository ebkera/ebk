"""
This script usually adds code that does not fit well into other categories
Some examples are 
	Error corrections
	
"""

def fix_ligand_names_TPDAC(ligand_name, complete_name = False):
	"""
	This is for the bringing the two naming conventions into the same standard.
		This hapened in TPDAc where there was a naming convention mismatch between the Sn atached ligands and the gas phase ligands.
	"""
	location_mapper_old_to_new_naming_convention = {"":"","A2":"A6", "A6": "A2", "A5":"A3","A3":"A5", "B2":"B3", "B3":"B2"}
	if complete_name:
		name = location_mapper_old_to_new_naming_convention[ligand_name]
		if name != "": return f"TPDAc_F_{name}"
		else: return f"TPDAc"
	return location_mapper_old_to_new_naming_convention[ligand_name]