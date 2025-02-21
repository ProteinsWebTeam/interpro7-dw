import oracledb
from datetime import datetime


def update_frozen_state(uri: str, is_frozen: bool | None, freeze_on: datetime | None):
    columns = []
    params = []
    if is_frozen is not None:
        columns.append("SET ACTIVE = :i")
        params.append("Y" if is_frozen else "N")

    if freeze_on is not None:
        columns.append("SET ACTIVE_FROM = :i")
        params.append(freeze_on)

    if columns:
        con = oracledb.connect(uri)
        cur = con.cursor()
        cur.execute(
            f"""
            UPDATE INTERPRO.PRONTO_STATES
            SET {','.join(columns)}
            WHERE NAME = 'FROZEN'
            """,
            params,
        )
        con.commit()
        cur.close()
        con.close()
        print("updated")
