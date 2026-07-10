import streamlit as st
from Bio.Seq import Seq

def analyze_dna(dna_string):
    """DNA 서열을 분석하는 함수"""
    dna_string = dna_string.upper().strip()
    
    valid_bases = set("ATGC")
    if not set(dna_string).issubset(valid_bases):
        return None
    
    seq = Seq(dna_string)
    
    count_a = dna_string.count("A")
    count_t = dna_string.count("T")
    count_g = dna_string.count("G")
    count_c = dna_string.count("C")
    
    gc_content = (count_g + count_c) / len(seq) * 100
    
    result = {
        "seq": seq,
        "length": len(seq),
        "count_a": count_a,
        "count_t": count_t,
        "count_g": count_g,
        "count_c": count_c,
        "gc_content": gc_content,
        "complement": seq.complement(),
        "reverse_complement": seq.reverse_complement(),
        "rna": seq.transcribe(),
    }
    
    try:
        result["protein"] = seq.translate()
    except Exception:
        result["protein"] = "번역 불가 (서열 길이 확인 필요)"
    
    return result


# ===== 웹 화면 구성 =====

st.title("🧬 DNA 서열 분석기")
st.write("DNA 서열을 입력하면 다양한 생물학적 정보를 분석해드려요.")

st.info("예시: ATGGCCATTGTAATG (A, T, G, C로만 입력)")

dna_input = st.text_input("DNA 서열을 입력하세요:")

if st.button("분석하기"):
    if dna_input:
        result = analyze_dna(dna_input)
        
        if result is None:
            st.error("⚠️ 올바른 DNA 서열이 아닙니다. A, T, G, C만 입력해주세요.")
        else:
            st.success(f"✅ 서열 길이: {result['length']} 염기")
            
            # 염기 구성 - 컬럼으로 예쁘게 배치
            st.write("### 염기 구성")
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("A (아데닌)", result["count_a"])
            col2.metric("T (티민)", result["count_t"])
            col3.metric("G (구아닌)", result["count_g"])
            col4.metric("C (사이토신)", result["count_c"])
            
            # GC 함량
            st.write("### GC 함량")
            st.progress(result["gc_content"] / 100)
            st.write(f"**{result['gc_content']:.2f}%** (GC 함량이 높을수록 DNA가 더 안정적입니다)")
            
            # 서열 정보들
            st.write("### 서열 분석 결과")
            st.code(f"입력 서열:     {result['seq']}\n"
                    f"상보적 서열:   {result['complement']}\n"
                    f"역상보 서열:   {result['reverse_complement']}\n"
                    f"RNA 전사:      {result['rna']}\n"
                    f"단백질 번역:   {result['protein']}")
    else:
        st.warning("DNA 서열을 입력해주세요!")