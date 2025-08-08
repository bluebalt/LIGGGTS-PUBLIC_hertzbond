# LIGGGHTS Input Script for Nanoindentation of C1 Anode Structure
units          micro
atom_style     granular

# 시뮬레이션 영역 정의 (150µm x 150µm x 76.5µm)
region         sample block 0 1.5e2  0 1.5e2  0 7.65e1 units box
create_box     1 sample  # 입자 타입 1개 (모든 입자 동일 재질)

# 입자 템플릿 정의 (5가지 크기 분포)
fix pts1 all particletemplate/sphere 15485863 atom_type 1 radius constant 1.5 density constant 2.2  # 지름 3µm
fix pts2 all particletemplate/sphere 15485867 atom_type 1 radius constant 2.5 density constant 2.2  # 지름 5µm
fix pts3 all particletemplate/sphere 32452843 atom_type 1 radius constant 3.5 density constant 2.2  # 지름 7µm
fix pts4 all particletemplate/sphere 32452867 atom_type 1 radius constant 4.5 density constant 2.2  # 지름 9µm
fix pts5 all particletemplate/sphere 49979687 atom_type 1 radius constant 6.0 density constant 2.2  # 지름 12µm

# 분포 정의: (지름 3:5:7:9:12 µm 입자의 수 비율 = 10:20:30:20:20)
fix pdd1 all particledistribution/discrete 67867967 5 pts1 0.10 pts2 0.20 pts3 0.30 pts4 0.20 pts5 0.20

# 입자 삽입: 한번에 전체 삽입, 영역 내 목표 입자수 8251개 (실험 C1 기준)
fix ins all insert/pack seed 49979693 distributiontemplate pdd1 insert_every once  overlapcheck yes particles_in_region 8251 region sample

# 중력 및 적분 설정
# gravity        9.81 vector 0.0 0.0 -1.0    # 중력 가속도 (아래로)
fix integrator all nve/sphere            # 병진+회전 운동 방정식 적분
# fix gforce all gravity 1.0 vector 0 0 -1  # 중력 적용 (배율 1.0 => 9.81 m/s^2)

# Pair 스타일 및 물성치 설정 (Hertzian + Bond 모델)
pair_style     hertzianbond 1.2e1        # 커스텀 pair_style (binderRatio=0.12 예시 입력)
pair_coeff     * * 0.52e9 0.30 1.0e10 0.5 5.0e6 5.0e6   # E=0.52 GPa, ν=0.30, bondStiffness=1e10, ζ=0.5, σn=τt=5e6 Pa

neighbor       7 bin                 # 이웃 리스트 간격 (정밀하게)
neigh_modify   every 1 delay 0 check yes

# 경계 조건: X,Y 주기적 / Z 하부 고정벽 (집전체)
boundary       p p f
fix bottom all wall/gran model hertzianbond plane 0.0 0.0 1.0 0.0 # z=0 평면 벽 (법선 +z 방향)

# 초기 동적 진정화: 입자 삽입 후 약간의 시간 진행
timestep       1.0e-7                   # 시간 스텝 (초) - 작은 값 (Hertz 접촉 안정 위해)
run            20000                    # 20000 스텝 (약 2e-3 s) -> 입자 정착 및 결합 생성

unfix ins                             # 삽입 fix 해제 (입자 생성 완료)

# 상부 인덴터(mesh) 설정: 평평한 원판 (직경 100µm) - STL 파일 불러오기
fix indenter all mesh/surface file flatPunch100um.stl type 1 scale 1.0 unit si
fix ind_move all move/mesh mesh indenter linear 0.0 0.0 -0.00015 # z-방향 0.15 mm/s로 이동

# 인덴터 하강 단계: 최대 압입 깊이 7.65e-6 m까지
run            500000  # (약 0.005 s) 인덴터 하강 + 하중 가압

# 인덴터 상승(제거) 단계: 원위치 복귀
fix ind_return all move/mesh mesh indenter linear 0.0 0.0 0.00015  # 반대 방향 속도로 상승
run            500000  # 인덴터 언로드

# 결과 출력 (예: 인덴터 힘-변위 곡선)
# indenter에 작용한 총 힘은 fix mesh/surface의 출력으로 얻을 수 있음 (에너지,힘 등의 누적)
fix_print      all print 1000 "$(step) $(fx_indenter) $(f_indenter[3])" file indent_force.csv

