#!/usr/bin/env python3
"""
Crear sistema.gro y sistema.top desde archivos ACPYPE consolidados
"""

def create_sistema_gro():
    """Combinar aceto.gro y pva.gro en sistema.gro"""
    print("📦 Creando sistema.gro...")
    
    # Leer archivos .gro
    with open('aceto.gro', 'r') as f:
        aceto_lines = f.readlines()
    
    with open('pva.gro', 'r') as f:
        pva_lines = f.readlines()
    
    # Extraer átomos (líneas 3 en adelante, excluyendo última línea de dimensiones)
    aceto_atoms = aceto_lines[2:-1]  # Sin título, sin número, sin dimensiones
    pva_atoms = pva_lines[2:-1]
    
    total_atoms = len(aceto_atoms) + len(pva_atoms)
    
    # Crear sistema.gro
    with open('sistema.gro', 'w') as f:
        f.write("Sistema Acetaminofen + PVA Polymer\n")
        f.write(f"{total_atoms:5d}\n")
        
        # Agregar átomos acetaminofen (residuo 1, ACE)
        for i, line in enumerate(aceto_atoms):
            # Formato: resnum resname atomname atomnum x y z
            atom_name = line[10:15].strip()
            x = float(line[20:28])
            y = float(line[28:36]) 
            z = float(line[36:44])
            f.write(f"{1:5d}{'ACE':<4s}{atom_name:>5s}{i+1:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        
        # Agregar átomos PVA (residuo 2, PVA)
        for i, line in enumerate(pva_atoms):
            atom_name = line[10:15].strip()
            x = float(line[20:28])
            y = float(line[28:36])
            z = float(line[36:44])
            f.write(f"{2:5d}{'PVA':<4s}{atom_name:>5s}{i+21:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
        
        # Dimensiones de caja (usar máximo)
        f.write("   40.00000   20.00000   10.00000\n")
    
    print(f"✅ sistema.gro creado ({total_atoms} átomos)")

def create_sistema_top():
    """Crear sistema.top con estructura jerárquica correcta"""
    print("📋 Creando sistema.top...")
    
    # Leer atomtypes consolidados
    with open('atomtypes_consolidados.txt', 'r') as f:
        atomtypes_content = f.read()
    
    # Crear sistema.top
    with open('sistema.top', 'w') as f:
        f.write("; Sistema: Acetaminofen + PVA Polymer\n")
        f.write("; Topología consolidada para resolver conflictos ACPYPE/GROMACS\n\n")
        
        f.write("[ defaults ]\n")
        f.write("; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n")
        f.write("1               2               yes             0.5     0.8333\n\n")
        
        f.write("[ atomtypes ]\n")
        f.write("; Atomtypes consolidados\n")
        f.write(atomtypes_content)
        f.write("\n")
        
        f.write("; Include molecule definitions (sin [atomtypes])\n")
        f.write('#include "aceto_clean.itp"\n')
        f.write('#include "pva_polymer_clean.itp"\n\n')
        
        f.write("[ system ]\n")
        f.write("; Name\n")
        f.write("Sistema Acetaminofen + PVA Polymer\n\n")
        
        f.write("[ molecules ]\n")
        f.write("; Compound        #mols\n")
        f.write("acetaminofen        1\n")
        f.write("pva_polymer         1\n")
    
    print("✅ sistema.top creado con estructura jerárquica correcta")

def verify_system():
    """Verificar archivos creados"""
    print("\n🔍 Verificando sistema creado...")
    
    # Verificar sistema.gro
    with open('sistema.gro', 'r') as f:
        gro_lines = f.readlines()
    
    print(f"📄 sistema.gro:")
    print(f"  Título: {gro_lines[0].strip()}")
    print(f"  Átomos: {gro_lines[1].strip()}")
    print(f"  Primer átomo ACE: {gro_lines[2][:15]}")
    print(f"  Primer átomo PVA: {gro_lines[22][:15] if len(gro_lines) > 22 else 'No encontrado'}")
    
    # Verificar sistema.top
    with open('sistema.top', 'r') as f:
        top_content = f.read()
    
    atomtypes_count = top_content.count('\n') - top_content.replace('[ atomtypes ]', '').count('\n')
    includes_count = top_content.count('#include')
    
    print(f"📄 sistema.top:")
    print(f"  Contiene [defaults]: {'✅' if '[ defaults ]' in top_content else '❌'}")
    print(f"  Contiene [atomtypes]: {'✅' if '[ atomtypes ]' in top_content else '❌'}")
    print(f"  Includes: {includes_count}")
    print(f"  Contiene [molecules]: {'✅' if '[ molecules ]' in top_content else '❌'}")

if __name__ == "__main__":
    print("🔧 CREANDO SISTEMA ACPYPE/GROMACS")
    print("=" * 35)
    
    create_sistema_gro()
    create_sistema_top()
    verify_system()
    
    print("\n✅ SISTEMA CREADO EXITOSAMENTE")
    print("Archivos generados:")
    print("  • sistema.gro")
    print("  • sistema.top")
