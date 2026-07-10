import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

def search_compound(name):
    """PubChem에서 화합물 이름으로 검색하는 함수"""
    try:
        compounds = pcp.get_compounds(name, 'name')
        
        if not compounds:
            print(f"⚠️ '{name}'에 대한 검색 결과가 없습니다.")
            return None
        
        compound = compounds[0]  # 첫 번째 결과 사용
        return compound
    
    except Exception as e:
        print(f"검색 중 오류 발생: {e}")
        return None


def analyze_compound(compound):
    """PubChem에서 가져온 화합물 정보를 분석하는 함수"""
    
    print(f"\n화합물 이름: {compound.iupac_name}")
    print(f"분자식: {compound.molecular_formula}")
    print(f"SMILES: {compound.canonical_smiles}")
    print(f"PubChem CID: {compound.cid}")
    print("-" * 40)
    
    # RDKit으로 추가 분석
    mol = Chem.MolFromSmiles(compound.canonical_smiles)
    
    if mol:
        mol_weight = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        num_h_donors = Descriptors.NumHDonors(mol)
        num_h_acceptors = Descriptors.NumHAcceptors(mol)
        
        print(f"분자량: {mol_weight:.2f} g/mol")
        print(f"LogP: {logp:.2f}")
        print(f"수소 결합 공여체: {num_h_donors}")
        print(f"수소 결합 수용체: {num_h_acceptors}")
        print("-" * 40)
        
        # Lipinski's Rule 체크
        rules_passed = 0
        if mol_weight <= 500:
            rules_passed += 1
        if num_h_donors <= 5:
            rules_passed += 1
        if num_h_acceptors <= 10:
            rules_passed += 1
        if logp <= 5:
            rules_passed += 1
        
        print(f"💊 Lipinski's Rule: {rules_passed}/4 통과")
        
        # 구조 이미지 저장
        img = Draw.MolToImage(mol, size=(400, 400))
        img.save("검색된_분자구조.png")
        print("📷 구조 이미지가 '검색된_분자구조.png'로 저장되었습니다.")


# 프로그램 시작
print("=== PubChem 화합물 검색기 ===")
print("영어 약물 이름을 입력하면 자동으로 정보를 가져와요.")
print()
print("예시: aspirin, ibuprofen, caffeine, penicillin, morphine")
print()

while True:
    drug_name = input("약물 이름을 입력하세요 (종료하려면 'q'): ")
    
    if drug_name.lower() == 'q':
        print("프로그램을 종료합니다.")
        break
    
    print("검색 중...")
    compound = search_compound(drug_name)
    
    if compound:
        analyze_compound(compound)
    
    print()