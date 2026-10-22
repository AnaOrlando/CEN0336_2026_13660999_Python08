#O início desse código é muito parecido com o script do exercício 1, porque também é preciso receber o arquivo multi FASTA, fazer a leitura e identificação das sequências

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

         seqs[current_seq] = "" #aqui temos que atribuir uma string vazia, que depois receberá o codon - isso é diferente do código da questão 1

      else: #se não é o título, significa que é a sequência do gene, então tem converter para maiúsculo para padronziar
          if current_seq is not None: #contando apenas as sequências  com título da sequência válido
             seqs[current_seq] += line.upper() #aqui junta todas as linhas da sequência e deixa em maiúscula

#A partir desse ponto, vamos fazer as etapas exclusivas do exercício 2

#Mandando o resultado para o arquivo Python_08.codons-frame-1.nt, conforme solicitado na questão

saida_filename = "Python_08.codons-frame-1.nt"

#Dividindo a sequência em códons

with open(saida_filename, "w") as out_file:
     for seq_id, sequence in seqs.items():
         codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)] #aqui separa de 3 em 3

         #formatando a parte de cima da sequência e os espaços entre os codons

         header = f"{seq_id}-frame-1-codons"
         codons_line = " ".join(codons)

         #escrevendo no arquivo de saída

         out_file.write(f"{header}\n")
         out_file.write(f"{codons_line}\n")

#Imprimindo a mensagem de finalização

print("Python_08.codons-frame-1.nt gerado com sucesso")
