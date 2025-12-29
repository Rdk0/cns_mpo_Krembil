import csv

def calculate_bbb_score(aro_r, ha, mw, hba, hbd, tpsa, pka):
    """
    Calculate BBB (Blood-Brain Barrier) Score based on molecular properties.
    reference https://pubs.acs.org/doi/10.1021/acs.jmedchem.9b01220
    Gupta et al. J. Med. Chem. 2019, 62, 21, 9824–9836
    The argument values (read from .csv) were calcualated using ACD labs Percepta but the ref.
    uses another software package.
    This is a first draft.  The paper incudes a python script but this a re-write based on
   authors' Excel file formulas (see the paper's supporting inforamtion stack)
    
    Args:
        aro_r: Number of Aromatic Rings (B3)
        ha: Number of Heavy Atoms (B4)
        mw: Molecular Weight (B5)
        hba: Number of Hydrogen Bond Acceptors (B6)
        hbd: Number of Hydrogen Bond Donors (B7)
        tpsa: Topological Polar Surface Area (B9)
        pka: pKa value (B10)
        
    Returns:
        BBB Score (C11)
    """
    
    # C3: Aromatic Rings score
    if aro_r == 0:
        c3 = 0.336376
    elif aro_r == 1:
        c3 = 0.816016
    elif aro_r == 2:
        c3 = 1
    elif aro_r == 3:
        c3 = 0.691115
    elif aro_r == 4:
        c3 = 0.199399
    else:
        c3 = 0
    
    # C4: Heavy Atoms score
    if ha <= 5 or ha > 45:
        c4 = 0
    else:
        c4 = (1/0.624231) * (0.0000443 * ha**3 - 0.004556 * ha**2 + 0.12775 * ha - 0.463)
    
    # B8: MWHBN calculation
    b8 = mw**(-0.5) * (hba + hbd)
    
    # C8: MWHBN score
    if b8 <= 0.05 or b8 > 0.45:
        c8 = 0
    else:
        c8 = (1/0.72258) * (26.733 * b8**3 - 31.495 * b8**2 + 9.5202 * b8 - 0.1358)
    
    # C9: TPSA score
    if tpsa == 0 or tpsa > 120:
        c9 = 0
    else:
        c9 = (1/0.9598) * (-0.0067 * tpsa + 0.9598)
    
    # C10: pKa score
    if pka <= 3 or pka > 11:
        c10 = 0
    else:
        c10 = (1/0.597488) * (0.00045068 * pka**4 - 0.016331 * pka**3 + 
                               0.18618 * pka**2 - 0.71043 * pka + 0.8579)
    
    # C11: BBB SCORE
    bbb_score = c3 + c4 + c8 * 1.5 + c9 * 2 + c10 * 0.5
    
    return bbb_score


def process_csv(input_file, output_file):
    """
    Read CSV file, calculate BBB score for each row, and save results.
    
    Args:
        input_file: Path to input CSV file
        output_file: Path to output CSV file
    """
    
    try:
        # Read the input CSV
        with open(input_file, 'r', newline='', encoding='utf-8') as infile:
            reader = csv.DictReader(infile)
            
            # Get the fieldnames and add 'bbb_score'
            fieldnames = reader.fieldnames
            if fieldnames is None:
                print("Error: CSV file appears to be empty or malformed")
                return
            
            output_fieldnames = list(fieldnames) + ['bbb_score']
            
            # Store all rows with calculated scores
            rows_with_scores = []
            
            # Process each row
            for row_num, row in enumerate(reader, start=2):  # start=2 because row 1 is header
                try:
                    # Extract values and convert to appropriate types
                    aro_r = int(float(row['aro_r']))
                    ha = int(float(row['ha']))
                    mw = float(row['mw'])
                    hba = int(float(row['hba']))
                    hbd = int(float(row['hbd']))
                    tpsa = float(row['tpsa'])
                    pka = float(row['pka'])
                    
                    # Calculate BBB score
                    score = calculate_bbb_score(aro_r, ha, mw, hba, hbd, tpsa, pka)
                    
                    # Add score to row
                    row['bbb_score'] = f"{score:.6f}"
                    rows_with_scores.append(row)
                    
                    print(f"Row {row_num}: BBB Score = {score:.6f}")
                    
                except (ValueError, KeyError) as e:
                    print(f"Warning: Skipping row {row_num} due to error: {e}")
                    continue
        
        # Write the output CSV
        with open(output_file, 'w', newline='', encoding='utf-8') as outfile:
            writer = csv.DictWriter(outfile, fieldnames=output_fieldnames)
            writer.writeheader()
            writer.writerows(rows_with_scores)
        
        print(f"\nSuccessfully processed {len(rows_with_scores)} rows")
        print(f"Results saved to: {output_file}")
        
    except FileNotFoundError:
        print(f"Error: Input file '{input_file}' not found")
    except Exception as e:
        print(f"Error processing CSV: {e}")


if __name__ == "__main__":
    # Specify input and output file paths
    input_csv = "\\folder_names\\acd.csv"
    output_csv = "\\folder_names\\acd_with_kremblin_bbb_scores.csv"
    
    print("BBB Score Calculator - CSV Processor")
    print("=" * 50)
    
    # Process the CSV
    process_csv(input_csv, output_csv)




    