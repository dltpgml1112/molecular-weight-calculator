import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
import io

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
        result["protein"] = "번역 불가"
    
    return result


# ===== 웹 화면 구성 =====

st.title("🧬📁 DNA 파일 분석기")
st.write("DNA 서열을 직접 입력하거나, FASTA 파일을 업로드해서 분석할 수 있어요.")

# 입력 방식 선택
input_method = st.radio("입력 방식을 선택하세요:", ["직접 입력", "파일 업로드"])

dna_sequence = None
sequence_name = None

if input_method == "직접 입력":
    dna_input = st.text_input("DNA 서열을 입력하세요:")
    if dna_input:
        dna_sequence = dna_input
        sequence_name = "직접 입력한 서열"

else:  # 파일 업로드
    uploaded_file = st.file_uploader("FASTA 파일을 업로드하세요 (.fasta, .fa, .txt)", 
                                       type=["fasta", "fa", "txt"])
    
    if uploaded_file is not None:
        # 파일 내용 읽기
        stringio = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        
        try:
            # FASTA 형식으로 읽기 시도
            records = list(SeqIO.parse(stringio, "fasta"))
            
            if records:
                # 여러 서열이 있으면 선택할 수 있게
                record_names = [r.id for r in records]
                selected = st.selectbox("분석할 서열을 선택하세요:", record_names)
                
                selected_record = next(r for r in records if r.id == selected)
                dna_sequence = str(selected_record.seq)
                sequence_name = selected_record.id
                
                st.success(f"'{sequence_name}' 서열을 불러왔습니다. (길이: {len(dna_sequence)} 염기)")
            else:
                st.error("FASTA 형식을 인식하지 못했습니다. 파일 형식을 확인해주세요.")
        except Exception as e:
            st.error(f"파일을 읽는 중 오류 발생: {e}")


# 분석하기 버튼
if dna_sequence:
    if st.button("분석하기"):
        result = analyze_dna(dna_sequence)
        
        if result is None:
            st.error("⚠️ 올바른 DNA 서열이 아닙니다. A, T, G, C만 포함되어야 합니다.")
        else:
            st.success(f"✅ [{sequence_name}] 서열 길이: {result['length']} 염기")
            
            st.write("### 염기 구성")
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("A", result["count_a"])
            col2.metric("T", result["count_t"])
            col3.metric("G", result["count_g"])
            col4.metric("C", result["count_c"])
            
            st.write("### GC 함량")
            st.progress(result["gc_content"] / 100)
            st.write(f"**{result['gc_content']:.2f}%**")
            
            st.write("### 서열 분석 결과")
            st.code(f"상보적 서열:   {result['complement']}\n"
                    f"역상보 서열:   {result['reverse_complement']}\n"
                    f"RNA 전사:      {result['rna']}\n"
                    f"단백질 번역:   {result['protein']}")