import sys

#Recebendo um arquivo multi-FASTA pela linha de comando

if len(sys.argv) > 1:
	filename: sys.arvg[1]
else:
	filename = input("Informe o nome do arquivo multi-FASTA:")

#Estruturando os dados para seguir a forma pedida: seqs[geneName][nucleotide]=count

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

         #Fazendo um dicionário para os nucleotídeos, começando com zero, para depois contar e somar a ele

         seqs[current_seq] = {"A" : 0, "T" : 0, "G" : 0, "C" : 0}

      else: #se não é o título, significa que é a sequência do gene, então tem converter para maiúsculo para padronziar
          if current_seq is not None: #contando apenas os nucleotídeos com título da sequência válido
             sequence = line.upper()

         #Contando os nucleotídeos de cada linha e somando ao dicionário

             for nucleotide in sequence:
                 if nucleotide in seqs[current_seq]:
                    seqs[current_seq][nucleotide] += 1

#Imprimindo o resultado. Coloquei também o nome das colunas, para facilitar interpretração

print("gene\t\tA\tT\tG\tC") #aqui foi necessário dar mais de um espaço entre o nome da sequência e os nucleotídeos devido ao espaçamento
for seq_name, counts in seqs.items():
    print(f"{seq_name}\t{counts['A']}\t{counts['T']}\t{counts['G']}\t{counts['C']}")
