import os
import sys
from Bio import SeqIO
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
import pysam

# 이전 phredqc에서 사용한 fastq 파일 ( 요루바족 여성 NA18489 )의 20번 염색체 엑솜 정렬파일을 사용.

# pysam : SAMTools C API의 파이썬 구현인 pysam 라이브러리를 활용.

# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/exome_alignment/NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam
# ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA18489/exome_alignment/NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam.bai
# reference : https://github.com/PacktPublishing/Bioinformatics-with-Python-Cookbook-Second-Edition/blob/master/Datasets.ipynb

# 참고자료 : http://www.incodom.kr/파이썬/라이브러리/pySam

# 코드를 작성하기전에 BAM파일을 확인한다. 언제나 분석전에 앞서 파일의 머리글(header)와 처음 몇 레코드를 살펴보는 습관을 기르자

bam = pysam.AlignmentFile(
    'NA18489.chrom20.ILLUMINA.bwa.YRI.exome.20121211.bam', 'rb')

headers = bam.header
'''
for record_type, records in headers.items():
    print (record_type)
    for i, record in enumerate(records):
        if type(record) == dict:
            print('\t%d' % (i + 1))
            for field, value in record.items():
                print('\t\t%s\t%s' % (field, value))
        else:
            print('\t\t%s' % record)
'''
# 위 헤더 정보는 텍스트로 정리했다.

for rec in bam:
    if rec.cigarstring.find('M') > -1 and \
            rec.cigarstring.find('S') > -1 and not \
            rec.is_unmapped and not rec.mate_is_unmapped:
        break
print(rec.query_name, rec.reference_id, bam.getrname(
    rec.reference_id), rec.reference_start, rec.reference_end)

print(rec.cigarstring)
print(rec.query_alignment_start, rec.query_alignment_end,
      rec.query_alignment_length)
print(rec.next_reference_id, rec.next_reference_start, rec.template_length)
print(rec.is_paired, rec.is_proper_pair, rec.is_unmapped, rec.mapping_quality)
print(rec.query_qualities)
print(rec.query_alignment_qualities)
print(rec.query_sequence)
