from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio import SeqIO

atp6_output_file = open("atp6_seq_features.fasta", "w")
cox1_output_file = open("cox1_seq_features.fasta", "w")
cytb_output_file = open("cytb_seq_features.fasta", "w")
atp8_output_file = open("atp8_seq_features.fasta", "w")
cox2_output_file = open("cox2_seq_features.fasta", "w")
cox3_output_file = open("cox3_seq_features.fasta", "w")
nd1_output_file = open("nd1_seq_features.fasta", "w")
nd2_output_file = open("nd2_seq_features.fasta", "w")
nd3_output_file = open("nd3_seq_features.fasta", "w")
nd4_output_file = open("nd4_seq_features.fasta", "w")
nd5_output_file = open("nd5_seq_features.fasta", "w")
nd6_output_file = open("nd6_seq_features.fasta", "w")
nd4l_output_file = open("nd4l_seq_features.fasta", "w")

for record in SeqIO.parse("genomic.gbff", "genbank"):
        for feature in record.features:
            if feature.type == 'gene' and 'pseudo' not in feature.qualifiers:
                if 'gene' in feature.qualifiers:
                    gene_name = feature.qualifiers['gene'][0]
                    if gene_name in ['ATP6', 'COX1', 'CYTB', 'ATP8', 'COX2', 'COX3', 'ND1', 'ND2', 'ND3', 'ND4', 'ND5', 'ND6', 'ND4L']:
                        location = feature.location
                        sequencias = record.seq[location.start:location.end]
                        if gene_name == 'ATP6':
                            atp6_output_file.write(f">{gene_name}\n")
                            atp6_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'COX1':
                            cox1_output_file.write(f">{gene_name}\n")
                            cox1_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'CYTB':
                            cytb_output_file.write(f">{gene_name}\n")
                            cytb_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'ATP8':
                            atp8_output_file.write(f">{gene_name}\n")
                            atp8_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'COX2':
                            cox2_output_file.write(f">{gene_name}\n")
                            cox2_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'COX3':
                            cox3_output_file.write(f">{gene_name}\n")
                            cox3_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'ND1':
                            nd1_output_file.write(f">{gene_name}\n")
                            nd1_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'ND2':
                            nd2_output_file.write(f">{gene_name}\n")
                            nd2_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'ND3':
                            nd3_output_file.write(f">{gene_name}\n")
                            nd3_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'ND4':
                            nd4_output_file.write(f">{gene_name}\n")
                            nd4_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'ND5':
                            nd5_output_file.write(f">{gene_name}\n")
                            nd5_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'ND6':
                            nd6_output_file.write(f">{gene_name}\n")
                            nd6_output_file.write(f"{sequencias}\n\n")
                        elif gene_name == 'ND4L':
                            nd4l_output_file.write(f">{gene_name}\n")
                            nd4l_output_file.write(f"{sequencias}\n\n")

# Fechar os arquivos de saída
atp6_output_file.close()
cox1_output_file.close()
cytb_output_file.close()
atp8_output_file.close()
cox2_output_file.close()
cox3_output_file.close()
nd1_output_file.close()
nd2_output_file.close()
nd3_output_file.close()
nd4_output_file.close()
nd5_output_file.close()
nd6_output_file.close()
nd4l_output_file.close()
atp6_output_file.close()
                    