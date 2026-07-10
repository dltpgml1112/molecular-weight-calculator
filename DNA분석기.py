from Bio.Seq import Seq

def analyze_dna(dna_string):
    """DNA 서열을 분석하는 함수"""
    
    # 대문자로 변환 (소문자로 입력해도 인식되게)
    dna_string = dna_string.upper()
    
    # 유효한 DNA인지 확인 (A, T, G, C만 있어야 함)
    valid_bases = set("ATGC")
    if not set(dna_string).issubset(valid_bases):
        print("⚠️ 올바른 DNA 서열이 아닙니다. A, T, G, C만 입력해주세요.")
        return
    
    # Biopython의 Seq 객체로 변환
    seq = Seq(dna_string)
    
    print(f"\n입력한 DNA 서열: {seq}")
    print(f"서열 길이: {len(seq)} 염기")
    print("-" * 40)
    
    # 각 염기 개수 세기
    count_a = dna_string.count("A")
    count_t = dna_string.count("T")
    count_g = dna_string.count("G")
    count_c = dna_string.count("C")
    
    print("염기 구성:")
    print(f"  A (아데닌): {count_a}개")
    print(f"  T (티민):   {count_t}개")
    print(f"  G (구아닌): {count_g}개")
    print(f"  C (사이토신): {count_c}개")
    print("-" * 40)
    
    # GC 함량 계산 (생물학적으로 중요한 지표)
    gc_content = (count_g + count_c) / len(seq) * 100
    print(f"GC 함량: {gc_content:.2f}%")
    print("(GC 함량이 높을수록 DNA가 더 안정적입니다)")
    print("-" * 40)
    
    # 상보적 서열 (Complementary strand)
    complement = seq.complement()
    print(f"상보적 서열:     {complement}")
    
    # 역상보 서열 (Reverse complement) - 실제 DNA 이중나선의 반대쪽 가닥
    reverse_complement = seq.reverse_complement()
    print(f"역상보 서열:     {reverse_complement}")
    print("-" * 40)
    
    # DNA를 RNA로 전사 (transcription)
    rna = seq.transcribe()
    print(f"RNA로 전사:      {rna}")
    
    # RNA를 단백질로 번역 (translation) - 3개씩 끊어서 아미노산으로 변환
    try:
        protein = seq.translate()
        print(f"단백질로 번역:   {protein}")
    except Exception as e:
        print("단백질 번역 중 오류 (서열 길이가 3의 배수가 아닐 수 있음)")


# 프로그램 시작
print("=== DNA 서열 분석기 ===")
print("A, T, G, C로 이루어진 DNA 서열을 입력하세요.")
print("예시: ATGGCCATTGTAATG")
print()

while True:
    user_input = input("DNA 서열을 입력하세요 (종료하려면 'q' 입력): ")
    
    if user_input.lower() == 'q':
        print("프로그램을 종료합니다.")
        break
    
    analyze_dna(user_input)
    print()