#!/usr/bin/env python

import os
import re
import subprocess
import sys


def align_fasta(fasta1_file, fasta2_file):
    """
    funkcja do alignowania sekwencji.  Sekwencje albo sa w oddzielnych plikach i wtedy laczymy je w jednym
    tmp file (usuwanym po wykonaniu funkcji), albo w jednym pliku. Do alignowania uzywamy maffta w strategii
    auto
    :param fasta1_file: str, sciezka do pliku fasta zawierajacego jedna lub wiecej sekwencji
    :param args: tupple, sciezka do dowolnej liczy plikow fasta
    :return: dict, slownik ktora jako klucze zawiera identyfikatory sekwencji z podanych pliku/plikow a jako klucze
    alignowano sekwencje
    """

    # Najpierw normalizujemy sciezki  isprawdzamy czy podane plik/pliki istnieja
    fasta1_file = [os.path.abspath(fasta1_file), os.path.abspath(fasta2_file)]
    stan1 = os.path.isfile(fasta1_file[0])
    stan2 = os.path.isfile(fasta1_file[1])
    if (not stan1) or (not stan2):
        raise Exception('Nie podano poprawnej sciezki do pliku!')

    for plik in fasta1_file:
        with open('tmp.fasta', 'a+') as f, open(plik, 'r') as f1:
            for line in f1:
                f.write(line.split()[0])
                f.write('\n')

    polecenie = ('mafft --auto --quiet --inputorder tmp.fasta')
    alignment = subprocess.Popen(polecenie, shell=True, stdout=subprocess.PIPE)
    # generowanie slownika
    slownik_alignmentu = {}
    for line in alignment.stdout:
        line = line.decode('UTF-8').rstrip()
        if len(re.findall(">", line)) > 0:
            klucz = line[:]
            if klucz not in slownik_alignmentu.keys():
                slownik_alignmentu[klucz] = ''
        else:
            slownik_alignmentu[klucz] = slownik_alignmentu[klucz] + line
    # pozbywamy sie tepowego pliku
    # os.remove('tmp.fasta')

    return slownik_alignmentu


if __name__ == '__main__':
    plik1 = sys.argv[1]  # to jest fasta z consensusem traktowana jako REFERENCJA
    plik2 = sys.argv[2]  # to jest fasta z dodanymi SV-s z manty traktowana jako TARGET
    out_name = sys.argv[3]  # output

    ref_name = open(plik1).readlines()[0].split()[0]
    target_name = open(plik2).readlines()[0].split()[0]

    sekwencja_with_n = ''
    with open(out_name, 'w') as f:
        aln = align_fasta(plik1, plik2)
        # teraz lecimy po
        target = aln[target_name]
        ref = aln[ref_name]
        for i in range(len(ref)):
            if target[i] == '-' and ref[i].lower() == 'n':
                sekwencja_with_n += '-'
            elif target[i] == '-' and ref[i].lower() != 'n':
                sekwencja_with_n += ref[i].upper()
            else:
                sekwencja_with_n += ref[i].upper()
        sekwencja_with_n = ''.join([element for element in sekwencja_with_n if element != '-'])
        f.write(ref_name + "_SV" + "\n")
        f.write(sekwencja_with_n)
        f.write('\n')
