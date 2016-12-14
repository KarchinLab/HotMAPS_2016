/* Adds the biomolecule num column to the output */
SELECT * FROM ( 
    SELECT DISTINCT SUBSTRING(gp.PDBId, 1, CHAR_LENGTH(gp.PDBId)-2) as PDBId, 
                    SUBSTRING(gp.PDBId, -1, 1) AS chain, 
                    pi.hugo,
                    bm.biomoleculeNo
    FROM PDB_Info as pi, Genome2PDB as gp, biomolecules as bm
    WHERE pi.pdbId=gp.PDBId and pi.modbase_filtered=1 and bm.pdbId=SUBSTRING(gp.PDBId, 1, CHAR_LENGTH(gp.PDBId)-2)
    UNION
    SELECT DISTINCT SUBSTRING(gp.PDBId, 1, CHAR_LENGTH(gp.PDBId)-2) as PDBId, 
                    SUBSTRING(gp.PDBId, -1, 1) AS chain, 
                    pi.hugo,
                    0 as biomoleculeNo
    FROM PDB_Info as pi, Genome2PDB as gp
    WHERE pi.pdbId=gp.PDBId and pi.modbase_filtered=1 and (pi.pdbId LIKE 'NP_%' or pi.pdbId LIKE 'ENSP%')
) o ORDER BY o.PDBId, o.hugo;
