from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def compare_compounds(smiles1, smiles2):
    """두 화합물의 구조적 유사도를 비교하는 함수"""
    
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)
    
    if mol1 is None or mol2 is None:
        print("⚠️ 올바른 SMILES 표기법이 아닙니다.")
        return
    
    # 분자 지문(fingerprint) 생성
    fp1 = AllChem.GetMorganFingerprintAsBitVect(mol1, radius=2, nBits=2048)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(mol2, radius=2, nBits=2048)
    
    # Tanimoto 유사도 계산
    similarity = DataStructs.TanimotoSimilarity(fp1, fp2)
    
    print(f"\n화합물 1: {smiles1}")
    print(f"화합물 2: {smiles2}")
    print("-" * 40)
    print(f"구조적 유사도 (Tanimoto): {similarity:.3f}")
    print(f"백분율로 표현: {similarity * 100:.1f}%")
    print("-" * 40)
    
    if similarity >= 0.7:
        print("🟢 매우 유사한 구조 (같은 계열의 약물일 가능성)")
    elif similarity >= 0.4:
        print("🟡 어느 정도 유사한 구조")
    else:
        print("🔴 구조적으로 많이 다름")


# 프로그램 시작
print("=== 화합물 유사도 비교기 ===")
print("두 개의 SMILES를 입력하면 구조적 유사도를 계산해드려요.")
print()
print("예시 화합물들:")
print("  아스피린:      CC(=O)OC1=CC=CC=C1C(=O)O")
print("  이부프로펜:    CC(C)CC1=CC=C(C=C1)C(C)C(=O)O")
print("  아세트아미노펜: CC(=O)NC1=CC=C(C=C1)O")
print("  카페인:        CN1C=NC2=C1C(=O)N(C(=O)N2C)C")
print()

while True:
    smiles1 = input("첫 번째 SMILES (종료하려면 'q'): ")
    if smiles1.lower() == 'q':
        print("프로그램을 종료합니다.")
        break
    
    smiles2 = input("두 번째 SMILES: ")
    if smiles2.lower() == 'q':
        print("프로그램을 종료합니다.")
        break
    
    compare_compounds(smiles1, smiles2)
    print()