from rdkit import Chem
from rdkit.Chem import Draw
import os
import pandas as pd
import argparse

def pareto_ligands_to_smiles_and_images(pdbid, experiment, epoch, num_gen, known_binding_site, output_dir=None):
    """
    Generate SMILES and images for ligands in the Pareto frontier.
    
    Args:
        pdbid: PDB ID (e.g., "4af3")
        experiment: Experiment name (e.g., "bmB")
        epoch: Epoch number
        num_gen: Number of generated molecules
        known_binding_site: Known binding site flag
        output_dir: Output directory (optional, will be created if not provided)
    """
    # Construct paths
    base_dir = f"{pdbid}/experiment_{experiment}_{epoch}_{num_gen}_{known_binding_site}_{pdbid}"
    pareto_csv = f"{base_dir}/pareto_front.csv"
    ligands_dir = f"{base_dir}/ligands"
    
    if output_dir is None:
        output_dir = f"{base_dir}/pareto_images"
    
    # Check if files exist
    if not os.path.exists(pareto_csv):
        print(f"Error: Pareto frontier file not found: {pareto_csv}")
        return
    
    if not os.path.exists(ligands_dir):
        print(f"Error: Ligands directory not found: {ligands_dir}")
        return
    
    # Create output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Read Pareto frontier CSV
    try:
        pareto_df = pd.read_csv(pareto_csv)
        print(f"Found {len(pareto_df)} molecules in Pareto frontier")
    except Exception as e:
        print(f"Error reading Pareto CSV: {e}")
        return
    
    # Get the ligand filenames from the Pareto frontier
    if 'filename' in pareto_df.columns:
        pareto_filenames = set(pareto_df['filename'].tolist())
    elif 'ligand' in pareto_df.columns:
        pareto_filenames = set(pareto_df['ligand'].tolist())
    else:
        print(f"Error: Could not find filename/ligand column in {pareto_csv}")
        print(f"Available columns: {pareto_df.columns.tolist()}")
        return
    
    smiles_list = []
    processed_count = 0
    
    # Process only ligands in Pareto frontier
    for filename in pareto_filenames:
        sdf_file = filename if filename.endswith('.sdf') else f"{filename}.sdf"
        sdf_path = os.path.join(ligands_dir, sdf_file)
        
        if not os.path.exists(sdf_path):
            print(f"Warning: SDF file not found: {sdf_path}")
            continue
        
        try:
            suppl = Chem.SDMolSupplier(sdf_path)
            for idx, mol in enumerate(suppl):
                if mol is None:
                    continue
                
                # Compute 2D coordinates for better depiction
                Chem.rdDepictor.Compute2DCoords(mol)
                smiles = Chem.MolToSmiles(mol)
                smiles_list.append((filename, smiles))
                
                print(f"Processing: {filename}")
                print(f"SMILES: {smiles}")
                
                # Use MolDraw2DCairo for higher quality images
                drawer = Draw.MolDraw2DCairo(500, 500)
                opts = drawer.drawOptions()
                opts.addAtomIndices = False
                opts.addStereoAnnotation = True
                opts.legendFontSize = 14
                opts.bondLineWidth = 2.0
                
                # Create legend with filename and SMILES (truncated if too long)
                smiles_display = smiles[:50] + "..." if len(smiles) > 50 else smiles
                legend = f"{os.path.splitext(filename)[0]}\n{smiles_display}"
                
                drawer.DrawMolecule(mol, legend=legend)
                drawer.FinishDrawing()
                img_bytes = drawer.GetDrawingText()
                
                img_path = os.path.join(output_dir, f"{os.path.splitext(filename)[0]}_mol_{idx+1}.png")
                with open(img_path, "wb") as img_file:
                    img_file.write(img_bytes)
                
                processed_count += 1
        except Exception as e:
            print(f"Error processing {sdf_file}: {e}")
    
    # Save SMILES to a text file
    smiles_file = os.path.join(output_dir, "pareto_smiles.txt")
    with open(smiles_file, "w") as f:
        f.write("filename\tsmiles\n")
        for filename, smiles in smiles_list:
            f.write(f"{filename}\t{smiles}\n")
    
    print(f"Processed {processed_count} Pareto frontier molecules. Images and SMILES saved to {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Generate SMILES and images for Pareto frontier ligands')
    parser.add_argument('--pdbid', type=str, required=True, help='PDB ID (e.g., 4af3)')
    parser.add_argument('--experiment', type=str, required=True, help='Experiment name (e.g., bmB)')
    parser.add_argument('--epoch', type=int, required=True, help='Epoch number')
    parser.add_argument('--num_gen', type=int, required=True, help='Number of generated molecules')
    parser.add_argument('--known_binding_site', type=str, required=True, help='Known binding site flag')
    parser.add_argument('--output_dir', type=str, required=False, help='Output directory (optional)')
    
    args = parser.parse_args()
    
    pareto_ligands_to_smiles_and_images(
        pdbid=args.pdbid,
        experiment=args.experiment,
        epoch=args.epoch,
        num_gen=args.num_gen,
        known_binding_site=args.known_binding_site,
        output_dir=args.output_dir
    )