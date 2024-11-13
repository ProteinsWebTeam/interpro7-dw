import oracledb


def get_interactions(url: str, raise_on_error: bool = True) -> dict:
    entries = {}
    try:
        con = oracledb.connect(url)
    except oracledb.OperationalError:
        if raise_on_error:
            raise
    else:
        cur = con.cursor()
        cur.execute(
            """
            SELECT 
              DISTINCT FX.PRIMARYID, C1.INTERACTION_AC, C1.INTERACTOR_AC, 
                UPPER(I1.SHORTLABEL), T1.SHORTLABEL, IX1.PRIMARYID, 
                C2.INTERACTOR_AC, UPPER(I2.SHORTLABEL), T2.SHORTLABEL, 
                IX2.PRIMARYID, PX.PRIMARYID
            FROM INTACT.IA_FEATURE_XREF FX
            INNER JOIN INTACT.IA_FEATURE F 
              ON FX.PARENT_AC = F.AC
            INNER JOIN INTACT.IA_RANGE R 
              ON F.AC = R.FEATURE_AC
            INNER JOIN INTACT.IA_COMPONENT C1 
              ON F.COMPONENT_AC = C1.AC
            INNER JOIN INTACT.IA_INTERACTOR I1 
              ON C1.INTERACTOR_AC = I1.AC
            INNER JOIN INTACT.IA_CONTROLLEDVOCAB T1 
              ON I1.INTERACTORTYPE_AC = T1.AC
            INNER JOIN INTACT.IA_INTERACTOR_XREF IX1 
              ON (I1.AC = IX1.PARENT_AC 
                AND IX1.DATABASE_AC='EBI-31'        -- uniprotkb
                AND IX1.QUALIFIER_AC='EBI-28')      -- identity
            INNER JOIN INTACT.IA_INT2EXP I2E 
              ON C1.INTERACTION_AC = I2E.INTERACTION_AC
            INNER JOIN INTACT.IA_EXPERIMENT E 
              ON I2E.EXPERIMENT_AC = E.AC
            INNER JOIN INTACT.IA_PUBLICATION P 
              ON (E.PUBLICATION_AC = P.AC 
                AND P.STATUS_AC='EBI-4326847')      -- released
            INNER JOIN INTACT.IA_PUBLICATION_XREF PX 
              ON (P.AC = PX.PARENT_AC 
                AND PX.DATABASE_AC='EBI-34' 
                  AND PX.QUALIFIER_AC='EBI-49940')  -- primary-reference
            INNER JOIN INTACT.IA_COMPONENT C2
              ON C1.INTERACTION_AC = C2.INTERACTION_AC 
                AND C1.AC < C2.AC
            INNER JOIN INTACT.IA_INTERACTOR I2 
              ON C2.INTERACTOR_AC = I2.AC
            INNER JOIN INTACT.IA_CONTROLLEDVOCAB T2 
              ON I2.INTERACTORTYPE_AC = T2.AC
            INNER JOIN INTACT.IA_INTERACTOR_XREF IX2 
              ON (I2.AC = IX2.PARENT_AC 
                AND IX2.DATABASE_AC='EBI-31'        -- uniprotkb
                AND IX2.QUALIFIER_AC='EBI-28')      -- identity
            WHERE FX.DATABASE_AC='EBI-65812'        -- interpro
            """
        )

        for row in cur:
            entry_acc = row[0]
            intact_id = row[1]
            interactor_1_id = row[2]
            protein_1_id = row[3]
            interactor_1_type = row[4]
            protein_1_acc = row[5]

            interactor_2_id = row[6]
            protein_2_id = row[7]
            interactor_2_type = row[8]
            protein_2_acc = row[9]

            try:
                pubmed_id = int(row[10])
            except (ValueError, TypeError):
                continue

            try:
                obj = entries[entry_acc]
            except KeyError:
                obj = entries[entry_acc] = []

            obj.append({
                "intact_id": intact_id,
                "pubmed_id": pubmed_id,
                "molecule_1": {
                    "accession": protein_1_acc,
                    "identifier": protein_1_id,
                    "type": interactor_1_type
                },
                "molecule_2": {
                    "accession": protein_2_acc,
                    "identifier": protein_2_id,
                    "type": interactor_2_type
                }
            })

        for k in entries:
            entries[k].sort(key=lambda x: x["intact_id"])

        cur.close()
        con.close()
    finally:
        return entries
