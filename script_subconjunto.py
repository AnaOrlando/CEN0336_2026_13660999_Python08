#O início (recebimento, leitura e processamento do arquivo FASTA) desse código é igual aos códigos dos exercícios anter>

import sys
import re

#Recebendo um arquivo multi-FASTA pela linha de comando

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    filename = input("Informe o nome do arquivo multi-FASTA: ").strip()

#Estruturando os dados

seqs = {}
current_seq = None

#Abrindo o arquivo para leitura e processamento

with open(filename, 'r') as file:
    for line in file:
        line = line.strip() #removendo quebra de linhas e espaços extras
        if not line:
            continue #pulando linhas vazias

        #Identificando o título da nova sequência (linha que começa com '>')

        if line.startswith(">"):
            current_seq = line[1:].split()[0] #aqui pega só o nome do gene (que vem depois do >) e exclui o >
            seqs[current_seq] = "" #aqui temos que atribuir uma string vazia, que depois receberá o codon - isso é dif>

        else:  #se não é o título, significa que é a sequência do gene, então tem converter para maiúsculo para padron>
            if current_seq is not None: #contando apenas as sequências  com título da sequência válido
                seqs[current_seq] += line.upper() #aqui junta todas as linhas da sequência e deixa em maiúscula

#Mandando o resultado para os arquivos solicitados

saida_filename_frames = "Python_08.codons-6frames.nt"
saida_filename_aa = "Python_08.translated.aa"
saida_filename_longest = "Python_08.translated-longest.aa"
saida_filename_subconjunto = "Python_08.orf-longest.nt"

#Criando uma tabela de complementariedade das bases

tabela_bases = str.maketrans("ATGC", "TACG")

#Criando uma tabela de tradução (dada no exercício, apenas copiada aqui no código)

tabela_traducao = {
    'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
    'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 'AGA':'R', 'AGG':'R',
    'AAT':'N', 'AAC':'N',
    'GAT':'D', 'GAC':'D',
    'TGT':'C', 'TGC':'C',
    'CAA':'Q', 'CAG':'Q',
    'GAA':'E', 'GAG':'E',
    'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G',
    'CAT':'H', 'CAC':'H',
    'ATT':'I', 'ATC':'I', 'ATA':'I',
    'TTA':'L', 'TTG':'L', 'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'AAA':'K', 'AAG':'K',
    'ATG':'M',
    'TTT':'F', 'TTC':'F',
    'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
    'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 'AGT':'S', 'AGC':'S',
    'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
    'TGG':'W',
    'TAT':'Y', 'TAC':'Y',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TAA':'*', 'TGA':'*', 'TAG':'*'
}

#Dividindo a sequência em códons

with open(saida_filename_frames, 'w') as out_file_frames, \
     open(saida_filename_aa, 'w') as out_file_aa, \
     open(saida_filename_longest, 'w') as out_file_longest, \
     open(saida_filename_subconjunto, 'w') as out_file_subconjunto:

    for seq_id, sequence in seqs.items():
        maior_peptideo = ""
        maior_frame = None
        maior_subconjunto = ""

        #Sequência direta (três primeiros frames de leitura, igual o exercício anterior)

        for frame in range(1, 4):
            start_index = frame - 1
            codons = [sequence[i:i+3] for i in range(start_index, len(sequence), 3)]

            #traduzindo apenas códons completos (de 3 bases)

            aminoacids_list = [tabela_traducao.get(c, 'X') for c in codons if len(c) == 3]

            #formatando a parte de cima da sequência e os espaços entre os códons

            header_frames = f">{seq_id}-frame-{frame}-codons"
            header_aa = f">{seq_id}-frame-{frame}-aminoacids"

            #escrevendo no arquivo de saída

            out_file_frames.write(f"{header_frames}\n{' '.join(codons)}\n")
            out_file_aa.write(f"{header_aa}\n{' '.join(aminoacids_list)}\n")

            #buscando peptídeos (M => Stop)

            aa_str = "".join(aminoacids_list)
            for match in re.finditer(r'M[^*]*\*', aa_str):
                idx_inicio = match.start()
                idx_fim = match.end()
                peptideo = aa_str[idx_inicio:idx_fim-1]

                if len(peptideo) > len(maior_peptideo):
                    maior_peptideo = peptideo
                    maior_frame = frame
                    maior_subconjunto = "".join(codons[idx_inicio:idx_fim])

        #Sequência complementar

        #Substituindo cada base pelo complemento e invertendo a string

        rev_sequence = sequence.translate(tabela_bases)[::-1]

        for frame in range(1, 4):
            frame_num = frame + 3
            start_index = frame - 1
            codons = [rev_sequence[i:i+3] for i in range(start_index, len(rev_sequence), 3)]

            #traduzindo apenas códons completos (de 3 bases)

            aminoacids_list = [tabela_traducao.get(c, 'X') for c in codons if len(c) == 3]

            #formatando a parte de cima da sequência e os espaços entre os códon

            header_frames = f">{seq_id}-frame-{frame_num}-codons"
            header_aa = f">{seq_id}-frame-{frame_num}-aminoacids"

            #escrevendo no arquivo de saída

            out_file_frames.write(f"{header_frames}\n{' '.join(codons)}\n")
            out_file_aa.write(f"{header_aa}\n{' '.join(aminoacids_list)}\n")

            #buscando por peptídeos (M => Stop) neste frame

            aa_str = "".join(aminoacids_list)
            for match in re.finditer(r'M[^*]*\*', aa_str):
                idx_inicio = match.start()
                idx_fim = match.end()
                peptideo = aa_str[idx_inicio:idx_fim-1]

                if len(peptideo) > len(maior_peptideo):
                    maior_peptideo = peptideo
                    maior_frame = frame_num
                    maior_subconjunto = "".join(codons[idx_inicio:idx_fim])

        #Salvando o maior peptídeo e sua sequência

        if maior_peptideo:
            header_longest = f">{seq_id}-frame-{maior_frame}-longest-peptide"
            out_file_longest.write(f"{header_longest}\n{maior_peptideo}\n")

            header_subconjunto = f">{seq_id}-frame-{maior_frame}-longest-subconjunto"
            out_file_subconjunto.write(f"{header_subconjunto}\n{maior_subconjunto}\n")

print("Python_08.codons-6frames.nt, Python_08.translated.aa, Python_08.translated-longest.aa e Python_08.orf-longest.nt gerados com sucesso!")
