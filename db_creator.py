"""
This is a table creator in your database for plasmid-search-engine
"""

import Addgene_parser
import psycopg2

DATABASE_CONFIG = {
    "database": "mydb",
    "user": "postgres",
    "password": "",
    "host": "127.0.0.1",  # local ip
    "port": 5432,
}


def get_connection():
    return psycopg2.connect(
        database=DATABASE_CONFIG.get('database'),
        user=DATABASE_CONFIG.get('user'),
        password=DATABASE_CONFIG.get('password'),
        host=DATABASE_CONFIG.get('host'),
        port=DATABASE_CONFIG.get('port'),
    )


def make_sql_style(plasmid: Addgene_parser.Plasmid) -> Addgene_parser.Plasmid:
    """Transform from Python's none to empty string and check for "'" sign
     and converts plasmid's sequnce into 16bit array"""
    for key, value in plasmid.__dict__.items():
        if isinstance(value, str):
            for letter in value:
                if letter == "'":
                    setattr(plasmid, f'{key}', "''".join(value.split("'")))  # removing "'" problem
    if isinstance(plasmid.sequence, (bytes, bytearray)):
        plasmid.sequence = plasmid.sequence.decode('utf-8')
    return plasmid


def create_table(id_list: list):
    print(f'Connecting to the database {DATABASE_CONFIG.get("database")}')
    with get_connection() as conn:
        Addgene_parser.PlasmidParser(id=id_list)
        for plasmid in Addgene_parser.PlasmidParser.plasmid_list:
            create_record(plasmid, conn)
            conn.commit()
        Addgene_parser.PlasmidParser.plasmid_list.clear()  # clear the cash of the parser

    print(f'Disconnecting from database {DATABASE_CONFIG.get("database")}')


def create_record(plasmid: Addgene_parser.Plasmid, conn):
    plasmid = make_sql_style(plasmid)

    cursor = conn.cursor()
    cursor.execute(
        "CREATE TABLE IF NOT EXISTS addgene_plasmids " +
        "(id INT PRIMARY KEY, name TEXT, size INT, backbone TEXT, vector_type TEXT, marker TEXT, resistance TEXT, " +
        "growth_t TEXT, growth_strain TEXT, growth_instructions TEXT, copy_num TEXT, gene_insert TEXT, url TEXT, "
        "sequence TEXT)")
    cursor.execute(
        f"INSERT INTO addgene_plasmids (id, name, size, backbone, vector_type, marker, resistance, growth_t, " +
        f"growth_strain, growth_instructions, copy_num, gene_insert, url, sequence)" +
        f"  VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)",
        (plasmid.id, plasmid.name, plasmid.size, plasmid.backbone, plasmid.vector_type, plasmid.marker,
         plasmid.resistance, plasmid.growth_t, plasmid.growth_strain, plasmid.growth_instructions,
         plasmid.copy_num, plasmid.gene_insert, plasmid.url, plasmid.sequence,))

    print(f"Plasmid {plasmid.id}, {plasmid.name} has been added to {DATABASE_CONFIG.get('database')}.")
    cursor.close()


def main():
    """Function that allows to parse all plasmids from Addgene"""
    current = 1
    step = 500
    start, end = current, (current + step)
    while end <= 150000:
        id_list = [i for i in range(start, end)]
        create_table(id_list)
        start += step
        end += step


if __name__ == '__main__':
    main()
