
#O início (recebimento, leitura e processamento do arquivo FASTA) desse código é igual aos códigos dos exercícios anteriores

import sys

#Recebendo um arquivo multi-FASTA pela linha de comando

if len(sys.argv) > 1:
        filename: sys.arvg[1]
else:
        filename = input("Informe o nome do arquivo multi-FASTA:")

#Estruturando os dados

seqs = {}
current_seq = None

#Abrindo o arquivo para leitura e processamento

with open(filename, 'r') as file:
     for line in file:
          line = line.strip() #removendo quebra de linhas e espaços extras
          if not line:
             continue #pulando linhas vazias

          #Identificando o título da nova sequência, ou seja, a linha que começa com '>', porque é a linha do nome do gene

          if line.startswith (">"):
             current_seq = line[1:].split()[0] #aqui pega só o nome do gene (que vem depois do >) e exclui o >

          #Fazendo um dicionário para os códons

             seqs[current_seq] = "" #aqui temos que atribuir uma string vazia, que depois receberá o codon - isso é dif>

          else: #se não é o título, significa que é a sequência do gene, então tem converter para maiúsculo para padron>
               if current_seq is not None: #contando apenas as sequências  com título da sequência válido
                  seqs[current_seq] += line.upper() #aqui junta todas as linhas da sequência e deixa em maiúscula

#Mandando o resultado para o arquivo Python_08.codons-6frames.nt, conforme solicitado na questão (aqui só muda o nome d>

saida_filename = "Python_08.codons-6frames.nt"

#Criando uma tabela de complementariedade das bases

tabela_bases = str.maketrans("ATGC", "TACG")

#Dividindo a sequência em códons

with open(saida_filename, 'w') as out_file:
    for seq_id, sequence in seqs.items():

        #Sequência direta (três primeiros frames de leitura, igual o exercício anterior)

        for frame in range(1, 4):
            start_index = frame - 1
            codons = [sequence[i:i+3] for i in range(start_index, len(sequence), 3)]

            #formatando a parte de cima da sequência e os espaços entre os códons

            header = f"{seq_id}-frame-{frame}-codons"
            codons_line = " ".join(codons)

            #escrevendo no arquivo de saída

            out_file.write(f"{header}\n")
            out_file.write(f"{codons_line}\n")

        #Sequência complementar

        #Substituindo cada base pelo complemento e invertendo a string

        rev_sequence = sequence.translate(tabela_bases)[::-1]
        for frame in range(1, 4):
            frame_num = frame + 3
            start_index = frame - 1
            codons = [rev_sequence[i:i+3] for i in range(start_index, len(rev_sequence), 3)]

            #formatando a parte de cima da sequência e os espaços entre os códons

            header = f"{seq_id}-frame-{frame_num}-codons"
            codons_line = " ".join(codons)

            #escrevendo no arquivo de saída

            out_file.write(f"{header}\n")
            out_file.write(f"{codons_line}\n")

print("Python_08.codons-6frames.nt gerado com sucesso")
