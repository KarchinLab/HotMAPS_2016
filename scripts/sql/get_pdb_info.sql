SELECT * FROM ( 
    SELECT DISTINCT SUBSTRING(gp.PDBId, 1, CHAR_LENGTH(gp.PDBId)-2) as PDBId, 
                    SUBSTRING(gp.PDBId, -1, 1) AS chain, 
                    pi.hugo 
    FROM PDB_Info as pi, Genome2PDB as gp 
    WHERE pi.pdbId=gp.PDBId and pi.modbase_filtered=1
) o ORDER BY o.PDBId, o.hugo;
