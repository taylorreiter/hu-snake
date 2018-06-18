#!/usr/bin/env python

import sys
from Bio.Blast import NCBIXML

OUT = open(sys.argv[1], 'w')
OUT.write("Query_Name\tQuery_Length\tAlignment_Title\tAlignment_ID\tAlignment_Def\tAlignment_Accession\tIdentitites\teValue\tBitscore\n")
for xml_file in sys.argv[2:]:
    result_handle = open(xml_file)
    blast_records = NCBIXML.parse(result_handle)
    for rec in blast_records:
        for alignment in rec.alignments:
            for hsp in alignment.hsps:
                fields = [rec.query, str(rec.query_length), alignment.title, alignment.hit_id,
                          alignment.hit_def, alignment.accession, str(hsp.identities), str(hsp.expect), str(hsp.bits)]
                OUT.write("\t".join(fields) + "\n")
OUT.close()
