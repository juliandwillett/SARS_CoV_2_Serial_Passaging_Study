#aligned to MN908947.3  (Wuhan ref - what Ensembl uses)
#Gene is the gene the position overlaps with, not considering
#  downstream/upstream gene variants

#VEP annotations done manually vs pulling from file to ensure picking highest impact 
#variants and it was strangely more time efficient.
#Verified on 29 Sep 2022

vep.df = data.frame(Pos=numeric(),Ref=character(),Alt=character(),
                    Gene=character(),AAPos=character(),
                    VEPAnnot=character(),VEPSeverity=character())
vep.df %<>% add_row(Pos=1,Ref='A',Alt='C',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=2,Ref='T',Alt='A',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=4,Ref='A',Alt='G',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=7,Ref='G',Alt='T',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifer') %>%
  add_row(Pos=7,Ref='G',Alt='C',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=8,Ref='G',Alt='C',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=10,Ref='T',Alt='G',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=174,Ref='G',Alt='T',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=210,Ref='G',Alt='T',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=240,Ref='TCG',Alt='TTG',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=241,Ref='C',Alt='T',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=314,Ref='AGTTT',Alt='AT',Gene='ORF1ab',AAPos='SL17-18M',VEPAnnot='Inframe Deletion',VEPSeverity='Moderate') %>%
  add_row(Pos=509,Ref='GGTCATGTTAT',Alt='GT',Gene='ORF1ab',AAPos='GHVM82-85V',VEPAnnot='Inframe Deletion',VEPSeverity='Moderate') %>%
  add_row(Pos=509,Ref='GGTCATGTTATGGTTG',Alt='GTGGTTG',Gene='ORF1ab',AAPos='GHVMVE82-87VVE',VEPAnnot='Inframe Deletion',VEPSeverity='Moderate') %>%
  add_row(Pos=543,Ref='A',Alt='G',Gene='ORF1ab',AAPos='E93G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=550,Ref='T',Alt='C',Gene='ORF1ab',AAPos='I95I',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=586,Ref='T',Alt='A',Gene='ORF1ab',AAPos='L107L',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=751,Ref='C',Alt='T',Gene='ORF1ab',AAPos='N162N',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=751,Ref='CACTAAACAT',Alt='TACTAAACAT',Gene='ORF1ab',AAPos='NTKH162-165NKTH',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=819,Ref='A',Alt='G',Gene='ORF1ab',AAPos='Y185C',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=919,Ref='A',Alt='G',Gene='ORF1ab',AAPos='Q218Q',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=929,Ref='A',Alt='G',Gene='ORF1ab',AAPos='I222V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=969,Ref='A',Alt='T',Gene='ORF1ab',AAPos='E235V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=1043,Ref='G',Alt='C',Gene='ORF1ab',AAPos='A260P',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=1059,Ref='C',Alt='T',Gene='ORF1ab',AAPos='T265I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=1131,Ref='A',Alt='G',Gene='ORF1ab',AAPos='E289G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=1310,Ref='TTGACT',Alt='TT',Gene='ORF1ab',AAPos='LT349-350X',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=1359,Ref='T',Alt='C',Gene='ORF1ab',AAPos='V365A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=1361,Ref='G',Alt='A',Gene='ORF1ab',AAPos='V366I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=1659,Ref='G',Alt='C',Gene='ORF1ab',AAPos='G465A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=1716,Ref='C',Alt='T',Gene='ORF1ab',AAPos='T484I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=1722,Ref='C',Alt='T',Gene='ORF1ab',AAPos='A486V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=2004,Ref='T',Alt='C',Gene='ORF1ab',AAPos='L580P',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=2005,Ref='C',Alt='T',Gene='ORF1ab',AAPos='L580L',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=2006,Ref='A',Alt='G',Gene='ORF1ab',AAPos='I581V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=2326,Ref='T',Alt='A',Gene='ORF1ab',AAPos='A687A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=2328,Ref='T',Alt='C',Gene='ORF1ab',AAPos='L688S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=2359,Ref='TAAA',Alt='GAAT',Gene='ORF1ab',AAPos='AK698-699AN',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=2362,Ref='A',Alt='T',Gene='ORF1ab',AAPos='K699N',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=2692,Ref='A',Alt='T',Gene='ORF1ab',AAPos='T809T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=2692,Ref='AAACAATAC',Alt='TAACAATAC',Gene='ORF1ab',AAPos='TNNT809-812TNNT',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=2786,Ref='A',Alt='C',Gene='ORF1ab',AAPos='I841L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=2791,Ref='T',Alt='C',Gene='ORF1ab',AAPos='T842T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=2911,Ref='T',Alt='C',Gene='ORF1ab',AAPos='T882T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=3024,Ref='TGTATTGTTCTTTCTACCCTCCAGATGAGGATGAAGAAGAAGG',Alt='TGTATTGTTCTTTTTACCCTCCAGATGAGGATGAAGAAGAAGG',Gene='ORF1ab',AAPos='MYCSFYPPDEDEEEG920-934MYCSFYPPDEDEEEG',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=3033,Ref='CTTTCT',Alt='CTTTTT',Gene='ORF1ab',AAPos='SFY923-925SFY',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=3037,Ref='C',Alt='T',Gene='ORF1ab',AAPos='F924F',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=3036,Ref='TCT',Alt='TTT',Gene='ORF1ab',AAPos='FY924-925FY',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=3042,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P926L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=3102,Ref='AATA',Alt='AA',Gene='ORF1ab',AAPos='QY946-947QX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=3230,Ref='G',Alt='T',Gene='ORF1ab',AAPos='G989C',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=3568,Ref='T',Alt='C',Gene='ORF1ab',AAPos='G1101G',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=3587,Ref='C',Alt='T',Gene='ORF1ab',AAPos='H1108Y',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=3669,Ref='ATGAAAATTTTA',Alt='AA',Gene='ORF1ab',AAPos='YENFN1135-1139*',VEPAnnot='Stop Gain/Frameshift',VEPSeverity='High') %>%
  add_row(Pos=3939,Ref='GAAAACAAG',Alt='GAAACAAG',Gene='ORF1ab',AAPos='RKQD1225-1228RNKX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=4115,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P1284S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=4213,Ref='T',Alt='C',Gene='ORF1ab',AAPos='A1316A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=4763,Ref='C',Alt='T',Gene='ORF1ab',AAPos='H1500Y',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=4861,Ref='T',Alt='C',Gene='ORF1ab',AAPos='D1532D',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=4981,Ref='A',Alt='G',Gene='ORF1ab',AAPos='T1572T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=5031,Ref='CATA',Alt='CA',Gene='ORF1ab',AAPos='TY1589-1590TX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=5157,Ref='CTTTTGA',Alt='CTTTGA',Gene='ORF1ab',AAPos='AFE1631-1633ALX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=5184,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P1640L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=5193,Ref='TGGGTA',Alt='TGGTA',Gene='ORF1ab',AAPos='LGR1643-1645LVX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=5224,Ref='TAAAAAGT',Alt='TAAAAATT',Gene='ORF1ab',AAPos='TKKW1653-1656TKNW',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=5353,Ref='A',Alt='G',Gene='ORF1ab',AAPos='Q1696Q',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=5500,Ref='A',Alt='G',Gene='ORF1ab',AAPos='K1745K',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=5673,Ref='CTTTTGTTA',Alt='CTTTGTTA',Gene='ORF1ab',AAPos='PFVM1803-1806PLLX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=5779,Ref='T',Alt='C',Gene='ORF1ab',AAPos='H1838H',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=5809,Ref='AGACGGTGCTTTACTTA',Alt='AA',Gene='ORF1ab',AAPos='DGALLT1849-1854T',VEPAnnot='Inframe Deletion',VEPSeverity='Moderate') %>%
  add_row(Pos=5909,Ref='TATAA',Alt='TA',Gene='ORF1ab',AAPos='YK1882-1883*',VEPAnnot='Stop Gained',VEPSeverity='High') %>%
  add_row(Pos=6070,Ref='C',Alt='T',Gene='ORF1ab',AAPos='I1935I',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=6207,Ref='A',Alt='C',Gene='ORF1ab',AAPos='K1981T',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=6297,Ref='GTTGTC',Alt='GTTTGTC',Gene='ORF1ab',AAPos='RCL2011-2013RLSX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=6449,Ref='C',Alt='T',Gene='ORF1ab',AAPos='L2062F',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=6515,Ref='T',Alt='C',Gene='ORF1ab',AAPos='L2084L',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=6813,Ref='C',Alt='T',Gene='ORF1ab',AAPos='T2183I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=6879,Ref='GTG',Alt='GCTG',Gene='ORF1ab',AAPos='SV2205-2206SCX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=6900,Ref='A',Alt='G',Gene='ORF1ab',AAPos='E2212G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=6989,Ref='T',Alt='C',Gene='ORF1ab',AAPos='S2242P',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=7119,Ref='C',Alt='T',Gene='ORF1ab',AAPos='S2285F',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=7224,Ref='CTTTTGGC',Alt='CTTTGGC',Gene='ORF1ab',AAPos='AFG2320-2322ALX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=7303,Ref='C',Alt='T',Gene='ORF1ab',AAPos='I2346I',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=7351,Ref='T',Alt='G',Gene='ORF1ab',AAPos='S2362S',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=7420,Ref='C',Alt='T',Gene='ORF1ab',AAPos='I2385I',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=7429,Ref='A',Alt='G',Gene='ORF1ab',AAPos='A2388A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=7471,Ref='CGGTT',Alt='CGTT',Gene='ORF1ab',AAPos='GC2403-2404VX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=7607,Ref='C',Alt='T',Gene='ORF1ab',AAPos='H2448Y',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=7623,Ref='T',Alt='C',Gene='ORF1ab',AAPos='V2453A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=8179,Ref='G',Alt='T',Gene='ORF1ab',AAPos='R2638R',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=8240,Ref='C',Alt='T',Gene='ORF1ab',AAPos='H2659Y',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=8431,Ref='T',Alt='C',Gene='ORF1ab',AAPos='S2722S',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=8443,Ref='A',Alt='C',Gene='ORF1ab',AAPos='R2726R',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=8466,Ref='A',Alt='T',Gene='ORF1ab',AAPos='K2734I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=8540,Ref='G',Alt='T',Gene='ORF1ab',AAPos='A2759S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=9059,Ref='T',Alt='C',Gene='ORF1ab',AAPos='Y2932H',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=9307,Ref='A',Alt='G',Gene='ORF1ab',AAPos='L3014L',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=9400,Ref='A',Alt='G',Gene='ORF1ab',AAPos='A3045A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=9621,Ref='A',Alt='G',Gene='ORF1ab',AAPos='D3119G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=9867,Ref='T',Alt='C',Gene='ORF1ab',AAPos='L3201P',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=9891,Ref='C',Alt='T',Gene='ORF1ab',AAPos='A3209V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=10320,Ref='TTAA',Alt='TTAG',Gene='ORF1ab',AAPos='LK3352-3353LR',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=10323,Ref='A',Alt='G',Gene='ORF1ab',AAPos='K3353R',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=10341,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P3359L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=10533,Ref='G',Alt='T',Gene='ORF1ab',AAPos='C3423F',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=10606,Ref='TT',Alt='CC',Gene='ORF1ab',AAPos='PF3447-3448PL',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=10809,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P3515L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=10834,Ref='C',Alt='T',Gene='ORF1ab',AAPos='A3523A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=10912,Ref='A',Alt='G',Gene='ORF1ab',AAPos='L3549L',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=11036,Ref='CTTTTAG',Alt='CTTTAG',Gene='ORF1ab',AAPos='LLV3591-3593L*X',VEPAnnot='Stop Gain/Frameshift',VEPSeverity='High') %>%
  add_row(Pos=11076,Ref='T',Alt='C',Gene='ORF1ab',AAPos='F3604S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=11101,Ref='ACCTTT',Alt='AT',Gene='ORF1ab',AAPos='PF3613-3614X',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=11104,Ref='T',Alt='C',Gene='ORF1ab',AAPos='P3613P',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=11134,Ref='T',Alt='C',Gene='ORF1ab',AAPos='A3623A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=11173,Ref='C',Alt='T',Gene='ORF1ab',AAPos='L3636L',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=11176,Ref='T',Alt='C',Gene='ORF1ab',AAPos='C3637C',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=11260,Ref='A',Alt='T',Gene='ORF1ab',AAPos='T3665T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=11275,Ref='TGATACTAGTTTGTCTGGTTTTA',Alt='TGATACTAGTTTGA',Gene='ORF1ab',AAPos='DTSLSGFK3671-3678DTSLK',VEPAnnot='Inframe Deletion',VEPSeverity='Moderate') %>%
  add_row(Pos=11275,Ref='TGATACTAGTTTGTCTGGTTTTA',Alt='TGTTTGA',Gene='ORF1ab',AAPos='DTSLSGFK3671-3678V*X',VEPAnnot='Stop Gained',VEPSeverity='High') %>%
  add_row(Pos=11287,Ref='GTCTGGTTTTA',Alt='GA',Gene='ORF1ab',AAPos='SGFK3675-3678K',VEPAnnot='Inframe Deletion',VEPSeverity='Moderate') %>%
  add_row(Pos=11303,Ref='A',Alt='G',Gene='ORF1ab',AAPos='K3680E',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=11355,Ref='C',Alt='T',Gene='ORF1ab',AAPos='A3697V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=11366,Ref='T',Alt='C',Gene='ORF1ab',AAPos='Y3701H',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=11418,Ref='T',Alt='C',Gene='ORF1ab',AAPos='V3718A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=11514,Ref='C',Alt='T',Gene='ORF1ab',AAPos='T3750I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=11776,Ref='C',Alt='T',Gene='ORF1ab',AAPos='A3837A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=12042,Ref='A',Alt='G',Gene='ORF1ab',AAPos='D3926G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=12075,Ref='A',Alt='G',Gene='ORF1ab',AAPos='N3937S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=12160,Ref='G',Alt='A',Gene='ORF1ab',AAPos='E3965E',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=12555,Ref='A',Alt='C',Gene='ORF1ab',AAPos='E4097A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=12636,Ref='G',Alt='T',Gene='ORF1ab',AAPos='W4124L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=12704,Ref='G',Alt='A',Gene='ORF1ab',AAPos='V4147I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=12775,Ref='T',Alt='C',Gene='ORF1ab',AAPos='A4170A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=12884,Ref='A',Alt='G',Gene='ORF1ab',AAPos='T4207A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=13145,Ref='T',Alt='C',Gene='ORF1ab',AAPos='C4294R',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=13168,Ref='C',Alt='T',Gene='ORF1ab',AAPos='H4301H',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=13192,Ref='A',Alt='G',Gene='ORF1ab',AAPos='T4309T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=13424,Ref='C',Alt='T',Gene='ORF1ab',AAPos='R4387C',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=13541,Ref='CTTTTGA',Alt='CTTTGA',Gene='ORF1ab',AAPos='AFD4426-4428ALX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=13830,Ref='T',Alt='C',Gene='ORF1ab',AAPos='A4522A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=14406,Ref='ACCTA',Alt='ACTTA',Gene='ORF1ab',AAPos='PPT4714-4716PLT',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=14406,Ref='ACCTACAA',Alt='ACTTACAA',Gene='ORF1ab',AAPos='PPTS4714-4717PLTS',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=14406,Ref='ACCTACAAGT',Alt='ACTTACAAGT',Gene='ORF1ab',AAPos='PPTS4714-4717PLTS',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=14407,Ref='CCTACAAGTTTTGGA',Alt='CTTACAAGTTTTGGA',Gene='ORF1ab',AAPos='PTSFG4715-4719LTSFG',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=14407,Ref='CC',Alt='CT',Gene='ORF1ab',AAPos='P4715L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=14408,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P4715L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=14496,Ref='T',Alt='C',Gene='ORF1ab',AAPos='G4744G',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=14574,Ref='T',Alt='C',Gene='ORF1ab',AAPos='P4770P',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=14776,Ref='G',Alt='A',Gene='ORF1ab',AAPos='G4838S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=14782,Ref='G',Alt='C',Gene='ORF1ab',AAPos='A4840P',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=15060,Ref='T',Alt='C',Gene='ORF1ab',AAPos='T4932T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=15199,Ref='GTAGTAATTG',Alt='GG',Gene='ORF1ab',AAPos='VVIG4979-4982GX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=15451,Ref='G',Alt='A',Gene='ORF1ab',AAPos='G5063S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=15460,Ref='T',Alt='C',Gene='ORF1ab',AAPos='Y5066H',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=15763,Ref='C',Alt='T',Gene='ORF1ab',AAPos='L5167L',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=15797,Ref='T',Alt='C',Gene='ORF1ab',AAPos='L5178P',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=15930,Ref='T',Alt='C',Gene='ORF1ab',AAPos='P5222P',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=15936,Ref='A',Alt='C',Gene='ORF1ab',AAPos='P5224P',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=16085,Ref='ATTTGTA',Alt='ATTGTA',Gene='ORF1ab',AAPos='HLY5274-5276HCX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=16275,Ref='A',Alt='T',Gene='ORF1ab',AAPos='S5337S',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=16289,Ref='C',Alt='T',Gene='ORF1ab',AAPos='A5342V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=16318,Ref='A',Alt='G',Gene='ORF1ab',AAPos='K5352E',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=16335,Ref='T',Alt='C',Gene='ORF1ab',AAPos='H5357H',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=16338,Ref='C',Alt='T',Gene='ORF1ab',AAPos='V5358V',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=16361,Ref='T',Alt='C',Gene='ORF1ab',AAPos='V5366A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=16466,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P5401L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=16630,Ref='CTTTTTGC',Alt='CTTTTTTGC',Gene='ORF1ab',AAPos='LFA5456-5458LFCX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=17821,Ref='C',Alt='A',Gene='ORF1ab',AAPos='P5853T',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=18471,Ref='T',Alt='C',Gene='ORF1ab',AAPos='D6069D',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=18586,Ref='T',Alt='C',Gene='ORF1ab',AAPos='F6108L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=19095,Ref='T',Alt='C',Gene='ORF1ab',AAPos='D6277D',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=19390,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P6376S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=19461,Ref='AACA',Alt='AA',Gene='ORF1ab',AAPos='T6400X',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=19683,Ref='A',Alt='C',Gene='ORF1ab',AAPos='E6473D',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=19884,Ref='CAAAAGAGAT',Alt='CAAGAGAT',Gene='ORF1ab',AAPos='KRD6541-6543KRX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=19905,Ref='T',Alt='C',Gene='ORF1ab',AAPos='H6547H',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=19983,Ref='CTTTTTTGA',Alt='CTTTTTGA',Gene='ORF1ab',AAPos='FFD6574-6576FLX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=19986,Ref='T',Alt='C',Gene='ORF1ab',AAPos='F6574F',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=20033,Ref='G',Alt='A',Gene='ORF1ab',AAPos='R6590H',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=20429,Ref='C',Alt='T',Gene='ORF1ab',AAPos='P6722L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=20594,Ref='C',Alt='T',Gene='ORF1ab',AAPos='T6777I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=20832,Ref='A',Alt='G',Gene='ORF1ab',AAPos='T6856T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=20937,Ref='G',Alt='A',Gene='ORF1ab',AAPos='T6891T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=21012,Ref='A',Alt='G',Gene='ORF1ab',AAPos='V6916V',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=21235,Ref='T',Alt='C',Gene='ORF1ab',AAPos='F6991L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=21255,Ref='G',Alt='T',Gene='ORF1ab',AAPos='A6997A',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=21557,Ref='C',Alt='T',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=21579,Ref='T',Alt='C',Gene='Spike',AAPos='V6A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=21611,Ref='AATCT',Alt='AATTT',Gene='Spike',AAPos='NL17-18NF',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=21612,Ref='ATC',Alt='ATT',Gene='Spike',AAPos='NL17-18NF',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=21614,Ref='C',Alt='T',Gene='Spike',AAPos='L18F',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=21618,Ref='C',Alt='G',Gene='Spike',AAPos='T19R',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=21651,Ref='A',Alt='C',Gene='Spike',AAPos='N30T',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=21796,Ref='GTTTGA',Alt='GTTTGC',Gene='Spike',AAPos='RFD78-80RFA',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=21801,Ref='A',Alt='C',Gene='Spike',AAPos='D80A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=22084,Ref='TTT',Alt='TCT',Gene='Spike',AAPos='PF174-175PL',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=22206,Ref='A',Alt='G',Gene='Spike',AAPos='D215G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=22227,Ref='C',Alt='T',Gene='Spike',AAPos='A222V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=22384,Ref='T',Alt='C',Gene='Spike',AAPos='T274T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=22433,Ref='T',Alt='C',Gene='Spike',AAPos='C291R',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=22459,Ref='A',Alt='G',Gene='Spike',AAPos='T299T',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=22480,Ref='C',Alt='T',Gene='Spike',AAPos='F306F',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=22595,Ref='A',Alt='G',Gene='Spike',AAPos='T345A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=22674,Ref='C',Alt='T',Gene='Spike',AAPos='S371F',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=22917,Ref='T',Alt='G',Gene='Spike',AAPos='L452R',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=22995,Ref='C',Alt='A',Gene='Spike',AAPos='T478K',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23011,Ref='TG',Alt='TA',Gene='Spike',AAPos='VE483-484VK',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23012,Ref='G',Alt='A',Gene='Spike',AAPos='E484K',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23062,Ref='TAA',Alt='TTA',Gene='Spike',AAPos='TN500-501TY',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23063,Ref='A',Alt='T',Gene='Spike',AAPos='N501Y',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23225,Ref='T',Alt='C',Gene='Spike',AAPos='S555P',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23246,Ref='TTCCAACAAT',Alt='TT',Gene='Spike',AAPos='FQQF562-565FX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=23277,Ref='C',Alt='T',Gene='Spike',AAPos='T572I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23291,Ref='C',Alt='T',Gene='Spike',AAPos='R577C',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23402,Ref='GAT',Alt='GGT',Gene='Spike',AAPos='D614G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23403,Ref='A',Alt='G',Gene='Spike',AAPos='D614G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23525,Ref='C',Alt='T',Gene='Spike',AAPos='H655Y',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23593,Ref='G',Alt='T',Gene='Spike',AAPos='Q677H',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23597,Ref='AATTCTCCTCG',Alt='ACTCGTCG',Gene='Spike',AAPos='NSPR679-682TRR',VEPAnnot='Protein altering variant',VEPSeverity='Moderate') %>%
  add_row(Pos=23597,Ref='AATTCTCC',Alt='AATTCTCG',Gene='Spike',AAPos='NSP679-681NSR',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23604,Ref='C',Alt='G',Gene='Spike',AAPos='P681R',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23606,Ref='C',Alt='T',Gene='Spike',AAPos='R682W',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23613,Ref='C',Alt='T',Gene='Spike',AAPos='A684V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23664,Ref='C',Alt='T',Gene='Spike',AAPos='A701V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=23664,Ref='CAGAAAATTC',Alt='TAGAAATTC',Gene='Spike',AAPos='AENS701-704VEIX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=23949,Ref='A',Alt='G',Gene='Spike',AAPos='D796G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=24237,Ref='C',Alt='T',Gene='Spike',AAPos='A892V',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=24375,Ref='T',Alt='C',Gene='Spike',AAPos='L938P',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25100,Ref='CAAAAAGAAATT',Alt='CAAAAAGGAATT',Gene='Spike',AAPos='QKEI1180-1183QKGI',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25107,Ref='A',Alt='G',Gene='Spike',AAPos='E1182G',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25333,Ref='T',Alt='C',Gene='Spike',AAPos='D1257D',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=25336,Ref='A',Alt='C',Gene='Spike',AAPos='E1258D',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25462,Ref='A',Alt='G',Gene='ORF3a',AAPos='T24A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25469,Ref='C',Alt='T',Gene='ORF3a',AAPos='S26L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25560,Ref='TCAG',Alt='TCAT',Gene='ORF3a',AAPos='FQ56-57FH',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25563,Ref='G',Alt='T',Gene='ORF3a',AAPos='Q57H',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25627,Ref='T',Alt='C',Gene='ORF3a',AAPos='F79L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25784,Ref='G',Alt='T',Gene='ORF3a',AAPos='W131L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=25792,Ref='C',Alt='A',Gene='ORF3a',AAPos='R134S',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=26031,Ref='G',Alt='A',Gene='ORF3a',AAPos='Q213Q',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=26160,Ref='TAATC',Alt='TATC',Gene='ORF3a',AAPos='NP257-258IX',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=26250,Ref='C',Alt='T',Gene='Envelope',AAPos='Y2Y',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=26456,Ref='C',Alt='T',Gene='Envelope',AAPos='P71L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=26607,Ref='C',Alt='T',Gene='Membrane',AAPos='L29F',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=26735,Ref='C',Alt='T',Gene='Membrane',AAPos='Y71Y',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=26767,Ref='T',Alt='C',Gene='Membrane',AAPos='I82T',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=26873,Ref='C',Alt='T',Gene='Membrane',AAPos='N117N',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=26927,Ref='A',Alt='G',Gene='Membrane',AAPos='E135E',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=26994,Ref='C',Alt='T',Gene='Membrane',AAPos='R158C',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=26995,Ref='G',Alt='A',Gene='Membrane',AAPos='R158H',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=27002,Ref='C',Alt='T',Gene='Membrane',AAPos='D160D',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=27044,Ref='A',Alt='G',Gene='Membrane',AAPos='R174R',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=27112,Ref='G',Alt='C',Gene='Membrane',AAPos='S197T',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=27179,Ref='G',Alt='A',Gene='Membrane',AAPos='L219L',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=27240,Ref='GATA',Alt='GA',Gene='ORF6',AAPos='I14X',VEPAnnot='Frameshift',VEPSeverity='High') %>%
  add_row(Pos=27393,Ref='C',Alt='T',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=27627,Ref='T',Alt='A',Gene='ORF7a',AAPos='R78R',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=27638,Ref='T',Alt='C',Gene='ORF7a',AAPos='V82A',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=27648,Ref='A',Alt='G',Gene='ORF7a',AAPos='K85K',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=27670,Ref='G',Alt='T',Gene='ORF7a',AAPos='V93F',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=27741,Ref='CAAAAGAAAGAC',Alt='CAAAAGAAAGAT',Gene='ORF7a',AAPos='LKRKT116-120LKRKI',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=27752,Ref='CA',Alt='TA',Gene='ORF7a',AAPos='T120I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=27752,Ref='CA',Alt='TC,TA',Gene='ORF7a',AAPos='T120I',VEPAnnot='Missense/Missense',VEPSeverity='Moderate/Moderate') %>%
  add_row(Pos=27752,Ref='CA',Alt='TC',Gene='ORF7a',AAPos='T120I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=27752,Ref='C',Alt='T',Gene='ORF7a',AAPos='T120I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=27925,Ref='C',Alt='T',Gene='ORF8',AAPos='T11I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=28068,Ref='G',Alt='T',Gene='ORF8',AAPos='E59*',VEPAnnot='Stop Gain',VEPSeverity='High') %>%
  add_row(Pos=28237,Ref='G',Alt='T',Gene='ORF8',AAPos='R115L',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=28242,Ref='GTTTTAGATTTCATCTAA',Alt='GTTTTAGATTTTATCTAA',Gene='ORF8',AAPos='VLDFI*117-122VLDFI*',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=28247,Ref='AGATTTCA',Alt='AA,AGATTTTA',Gene='ORF8',AAPos='DF119-121I/DFI119-121DFI',VEPAnnot='Inframe Deletion/Synonymous',VEPSeverity='Moderate/Low') %>%
  add_row(Pos=28247,Ref='AGATTTCA',Alt='AGATTTTA',Gene='ORF8',AAPos='DFI119-121DFI',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=28247,Ref='AGATTTCA',Alt='AA',Gene='ORF8',AAPos='DF119-121I',VEPAnnot='Inframe Deletion',VEPSeverity='Moderate') %>%
  add_row(Pos=28250,Ref='TTTCAT',Alt='TTTTAT',Gene='ORF8',AAPos='DFI119-121DFI',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=28251,Ref='TTCAT',Alt='TTTAT',Gene='ORF8',AAPos='FI120-121FI',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=28252,Ref='TCA',Alt='TTA',Gene='ORF8',AAPos='FI120-121FI',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=28253,Ref='CA',Alt='TA',Gene='ORF8',AAPos='FI120-121FI',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=28253,Ref='CA',Alt='TC',Gene='ORF8',AAPos='FI120-121FL',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=28253,Ref='CA',Alt='TC,TA',Gene='ORF8',AAPos='FI120-121FI/FI120-121FL',VEPAnnot='Synonymous/Missense',VEPSeverity='Low/Moderate') %>%
  add_row(Pos=28253,Ref='C',Alt='T',Gene='ORF8',AAPos='F120F',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=28270,Ref='TAAAATG',Alt='TAAATG',Gene='Nucleoprotein',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=28320,Ref='C',Alt='T',Gene='Nucleoprotein',AAPos='T16M',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=28368,Ref='G',Alt='A',Gene='Nucleoprotein',AAPos='R32H',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=28631,Ref='G',Alt='A',Gene='Nucleoprotein',AAPos='G120R',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=28881,Ref='G',Alt='T',Gene='Nucleoprotein',AAPos='R203M',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=28887,Ref='C',Alt='T',Gene='Nucleoprotein',AAPos='T205I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=29118,Ref='C',Alt='T',Gene='Nucleoprotein',AAPos='T282I',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=29200,Ref='C',Alt='T',Gene='Nucleoprotein',AAPos='P309P',VEPAnnot='Synonymous',VEPSeverity='Low') %>%
  add_row(Pos=29402,Ref='G',Alt='T',Gene='Nucleoprotein',AAPos='D377Y',VEPAnnot='Missense',VEPSeverity='Moderate') %>%
  add_row(Pos=29742,Ref='G',Alt='T',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=29746,Ref='AGTACGATCGAGTG',Alt='AGTATCGAGTG',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=29748,Ref='TACGA',Alt='TA',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=29748,Ref='TACGATCGAGTGTACAGTGAA',Alt='TATCGAGTGTACAGTGAA',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=29821,Ref='T',Alt='G',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=29868,Ref='G',Alt='T',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') %>%
  add_row(Pos=29868,Ref='G',Alt='A',Gene='Noncoding',AAPos='NA',VEPAnnot='Noncoding',VEPSeverity='Modifier') 
for (row in 1:nrow(vep.df)) {
  if (vep.df$VEPAnnot[[row]] == 'Frameshift') vep.df$VEPSeverity[[row]] = 'High'
  else if (vep.df$VEPAnnot[[row]] == 'Missense') vep.df$VEPSeverity[[row]] = 'Moderate'
  else if (vep.df$VEPAnnot[[row]] == 'Missense/Missense') vep.df$VEPSeverity[[row]] = 'Moderate'
  else if (vep.df$VEPAnnot[[row]] == 'Synonymous') vep.df$VEPSeverity[[row]] = 'Low'
  else if (vep.df$VEPAnnot[[row]] == 'Noncoding') vep.df$VEPSeverity[[row]] = 'Modifier'
}

getVEPData = function(df,print.vep.for.vep) { #print lines to plug into VEP or above data frame, saving time.
  out.df = df
  lines.for.vep = character()
  for (item in 1:nrow(df)) {
    curr.row = df[item,]
    m = which(vep.df$Pos==curr.row$Pos & vep.df$Ref==curr.row$Ref &
                vep.df$Alt==curr.row$Alt)
    lines.for.vep %<>% append(glue('MN908947.3  {curr.row$POS}  .  {curr.row$Ref}  {curr.row$Alt}'))
    if (length(m)==0 & print.vep.for.vep) { print(glue('MN908947.3  {curr.row$POS}  .  {curr.row$Ref}  {curr.row$Alt}')) ; next }
    else if (length(m)==0 & !print.vep.for.vep) { print(glue("add_row(Pos={curr.row$Pos},Ref='{curr.row$Ref}',Alt='{curr.row$Alt}',Gene='NA',AAPos='NA',VEPAnnot='NA',VEPSeverity='NA') %>%")) ; next }
    vep.match = vep.df[m,]
    out.df$Gene[[item]] = vep.match$Gene
    out.df$AAChange[[item]] = vep.match$AAPos
    out.df$VEPConsequence[[item]] = vep.match$VEPAnnot
    out.df$VEPSeverity[[item]] = vep.match$VEPSeverity
  }
  if (print.vep.for.vep) return(lines.for.vep)
  else return(out.df)
}
getDifferentVariantsList = function(df1,samples,return.justdf2) {
  out.df = samples[[1]][0,]
  for (item in 1:length(samples)) {
    tmp = getDifferentVariants(df1,samples[[item]],return.justdf2)
    out.df %<>% add_row(tmp)
    print(glue('Novel variants for sample {item}: {nrow(tmp)}'))
  }
  return(out.df %>% distinct(Pos,Ref,Alt,.keep_all = T) %>%
           dplyr::arrange(Pos))
}
getDifferentVariants = function(df1,df2,return.justdf2) {
  overlap.1 = which(df1$Pos %notin% df2$Pos)
  overlap.2 = which(df2$Pos %notin% df1$Pos)
  if (!return.justdf2) {
    df.out = df1[overlap.1,] %>% add_row(df2[overlap.2,])
    return(df.out %>% dplyr::arrange(Pos))
  }else{
    return(df2[overlap.2,] %>% dplyr::arrange(Pos))
  }
}
printVariantVariable = function(df) {
  plt.df = data.frame(Pos=numeric(),Passage=character())
  for (item in 1:length(df)) { #go by passage
    for (item2 in 1:length(df[[item]])) { #go by sample
      plt.df %<>% add_row(Pos=df[[item]][[item2]][['Pos']],
                          Passage=as.character(item))
    }
  }
  # return(plt.df)
  ggplot(plt.df,aes(x=Pos,fill=Passage)) + geom_histogram(position='identity') +
    theme_bw() + xlab('Position') + theme(text=element_text(size=20)) +
    ylab('')
  return(plt.df)
  # ggplot(plt.df,aes(x=Pos,y=0,color=Gene,shape=Gene,fill=Gene)) + scale_shape_manual(values=18:(18+length(unique(df$Gene)))) +
  #   geom_point(size=6) + xlim(c(0,30000)) +
  #   theme_bw() + theme(text=element_text(size=20),axis.text.y=element_blank(),axis.ticks.y=element_blank()) + 
  #   xlab('Position') + ylab('')
}
printVariantGeneInvolvement = function(df) { #take in 
  df$Gene = str_replace_all(df$Gene,'Noncoding','NC')
  df$Gene = str_replace_all(df$Gene,'Membrane','M')
  df$Gene = str_replace_all(df$Gene,'Nucleoprotein','NP')
  ggplot(df,aes(x=Gene,fill=Gene)) + geom_bar(aes(y=..count../sum(..count..)),show.legend = F) + 
    theme_bw() + theme(text=element_text(size=20)) + 
    ylab('Proportion') + ylim(c(0,0.5))
}
printVEPAnnotations = function(df) {
  df[which(is.na(df$VEPConsequence)),'VEPConsequence'] = 'Noncoding'
  ggplot(df,aes(y=VEPConsequence,fill=VEPConsequence)) + geom_bar(aes(x=..count../sum(..count..)),show.legend = F) + 
    theme_bw() + theme(text=element_text(size=28)) + xlab('Proportion') + 
    ylab('Functional Consequence') + xlim(c(0,0.8))
}
