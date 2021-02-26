import os
import gzip
from Bio import SeqIO
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
#%matplotlib inline

recs = SeqIO.parse(gzip.open('SRR003265.filt.fastq.gz', 'rt', encoding='utf-8'), 'fastq')
# fastq 파일을 파싱

#변수의 자료형을 정해준다
# 사전형 데이터에 기본값을 정의하고 만약 키값이 없더라도 오류를 출력하지
# 않고 미리 지정한 기본값을 출력한다.

cnt_qual = defaultdict(int)

# 레코드를 순회하면서 Phred 점수의 주석을 저장한다.
# 코드에서 사용한 fastq 파일은 이미 24개까지 필터링이 되었기 때문에
# 처음 24개는 무시하고 이후부터 저장한다.
'''
for rec in recs:
    for i, qual in enumerate(rec.letter_annotations['phred_quality']):
        if i < 25:
            continue
        cnt_qual[qual] += 1
tot = sum(cnt_qual.values())

for qual, cnt in cnt_qual.items():
    print('%d: %.2f %d' % (qual, 100. * cnt / tot, cnt))

'''
## 시각화

qual_pos = defaultdict(list)
for rec in recs:
    for i , qual in enumerate(rec.letter_annotations['phred_quality']):
        if i < 25 or qual == 40:
            continue
        pos = i + 1
        qual_pos[pos].append(qual)
        
vps = []
poses = list(qual_pos.keys())
poses.sort()

for pos in poses:
    vps.append(qual_pos[pos])
    
fig, ax = plt.subplots(figsize=(16,9))
sns.boxplot(data=vps, ax=ax)
ax.set_xticklabels(
    [str(x) for x in range(26, max(qual_pos.keys()) + 1)])

#plt.show()
plt.savefig('phred_quality.png')