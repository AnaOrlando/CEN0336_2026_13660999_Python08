#Grande parte desse código é igual ao script do exercício 2

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

             seqs[current_seq] = "" #aqui temos que atribuir uma string vazia, que depois receberá o codon - isso é diferen>

          else: #se não é o título, significa que é a sequência do gene, então tem converter para maiúsculo para padronziar
               if current_seq is not None: #contando apenas as sequências  com título da sequência válido
                  seqs[current_seq] += line.upper() #aqui junta todas as linhas da sequência e deixa em maiúscula

#Mandando o resultado para o arquivo Python_08.codons-3frames.nt, conforme solicitado na questão (aqui só muda o nome do arquivo, em relação ao script do exercício 2)

saida_filename = "Python_08.codons-3frames.nt"

#Dividindo os códons nos três primeiros quadros de leitura - esse trecho, antes de formatar a parte de cima e os espaços, é diferente do script do exercício 2

with open(saida_filename, 'w') as out_file:
     for seq_id, sequence in seqs.items():
         for frame in range(1, 4):
             start_index = frame - 1  #frame 1 = posição 0; frame 2 = posição 1 e frame 3 = posição 2
             codons = [sequence[i:i+3] for i in range(start_index, len(sequence), 3)] #aqui separa de 3 em 3, mas diferente do script do exercício 2, considera o index de início

             #formatando a parte de cima da sequência e os espaços entre os códons

             header = f"{seq_id}-frame-{frame}-codons"
             codons_line = " ".join(codons)

             #escrevendo no arquivo de saída

             out_file.write(f"{header}\n")
             out_file.write(f"{codons_line}\n")

#Imprimindo a mensagem de finalização

print("Python_08.codons-3frames.nt foi gerado com sucesso")
