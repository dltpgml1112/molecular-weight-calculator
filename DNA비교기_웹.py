import streamlit as st
from Bio.Seq import Seq

def compare_sequences(seq1, seq2):
    """두 DNA 서열을 비교하는 함수"""
    seq1 = seq1.upper().strip()
    seq2 = seq2.upper().strip()
    
    valid_bases = set("ATGC")
    if not set(seq1).issubset(valid_bases) or not set(seq2).issubset(valid_bases):
        return None
    
    # 짧은 쪽 길이에 맞춰서 비교
    min_length = min(len(seq1), len(seq2))
    
    matches = 0
    diff_positions = []
    
    for i in range(min_length):
        if seq1[i] == seq2[i]:
            matches += 1
        else:
            diff_positions.append(i + 1)  # 사람이 읽기 쉽게 1부터 시작
    
    similarity = (matches / min_length) * 100
    
    return {
        "seq1": seq1,
        "seq2": seq2,
        "min_length": min_length,
        "matches": matches,
        "similarity": similarity,
        "diff_positions": diff_positions,
        "diff_count": len(diff_positions)
    }


# ===== 웹 화면 구성 =====

st.title("🧬🔬 DNA 서열 비교기")
st.write("두 개의 DNA 서열을 비교해서 얼마나 비슷한지 계산해드려요.")

st.info("예시: 사람과 침팬지 유전자, 바이러스 변이주 등을 비교해보세요!")

col1, col2 = st.columns(2)

with col1:
    seq1_input = st.text_input("첫 번째 DNA 서열:", key="seq1")

with col2:
    seq2_input = st.text_input("두 번째 DNA 서열:", key="seq2")

if st.button("비교하기"):
    if seq1_input and seq2_input:
        result = compare_sequences(seq1_input, seq2_input)
        
        if result is None:
            st.error("⚠️ 올바른 DNA 서열이 아닙니다. A, T, G, C만 입력해주세요.")
        else:
            # 유사도 표시
            st.write("### 유사도 결과")
            st.metric("서열 유사도", f"{result['similarity']:.2f}%")
            st.progress(result['similarity'] / 100)
            
            st.write(f"비교한 길이: **{result['min_length']}** 염기 (짧은 쪽 기준)")
            st.write(f"일치하는 염기: **{result['matches']}** 개")
            st.write(f"다른 염기: **{result['diff_count']}** 개")
            
            # 다른 위치 표시
            if result['diff_count'] > 0:
                st.write("### 차이나는 위치 (돌연변이 후보)")
                positions_str = ", ".join(str(p) for p in result['diff_positions'][:20])
                st.code(positions_str)
                if result['diff_count'] > 20:
                    st.caption(f"... 외 {result['diff_count'] - 20}개 더 있음")
            else:
                st.success("두 서열이 완전히 동일합니다!")
            
            # 서열 나란히 비교해서 보여주기
            st.write("### 서열 비교")
            st.code(f"서열1: {result['seq1'][:result['min_length']]}\n"
                    f"서열2: {result['seq2'][:result['min_length']]}")
    else:
        st.warning("두 서열 모두 입력해주세요!")
    