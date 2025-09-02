#!/usr/bin/env python3
"""
Generar todos los archivos .mdp necesarios para simulación MD completa
Sistema: Acetaminofen + PVA Polymer
"""

def create_temporal_mdp():
    """MDP temporal para verificación de topología y genion"""
    content = """; Archivo MDP temporal para verificación de topología
; Solo para gmx grompp y gmx genion

integrator      = steep         ; Algoritmo (steep = minimización)
nsteps          = 0             ; Sin pasos (solo verificación)
emtol           = 1000.0        ; Tolerancia fuerza (kJ/mol/nm)
emstep          = 0.01          ; Tamaño paso inicial (nm)

; Parámetros de salida mínimos
nstlog          = 1             ; Frecuencia log
nstenergy       = 1             ; Frecuencia energía

; Control de temperatura - DESACTIVADO
tcoupl          = no            ; Sin termostato

; Control de presión - DESACTIVADO  
pcoupl          = no            ; Sin barostato

; Parámetros de fuerza (mínimos para verificación)
coulombtype     = Cut-off       ; Electrostática simple
rcoulomb        = 1.0           ; Radio corte coulomb (nm)
rvdw            = 1.0           ; Radio corte van der Waals (nm)
rlist           = 1.0           ; Radio lista vecinos (nm)

; Condiciones de frontera
pbc             = xyz           ; Condiciones frontera periódicas"""

    with open('temporal.mdp', 'w') as f:
        f.write(content)
    print("✅ temporal.mdp creado")

def create_ions_mdp():
    """MDP específico para genion (añadir iones)"""
    content = """; Archivo MDP para genion - Añadir iones
; Parámetros mínimos para genion

integrator      = steep         ; Algoritmo de minimización
nsteps          = 0             ; Sin pasos de simulación
emtol           = 1000.0        ; Tolerancia (kJ/mol/nm)
emstep          = 0.01          ; Tamaño paso (nm)

; Salida mínima
nstlog          = 1             ; Frecuencia escritura log
nstenergy       = 1             ; Frecuencia escritura energía

; Sin control temperatura/presión
tcoupl          = no            ; Sin termostato
pcoupl          = no            ; Sin barostato

; Parámetros electrostáticos (importantes para iones)
coulombtype     = PME           ; Particle Mesh Ewald
rcoulomb        = 1.0           ; Radio corte coulomb (nm)
rvdw            = 1.0           ; Radio corte vdW (nm)
rlist           = 1.0           ; Radio lista vecinos (nm)

; Condiciones frontera
pbc             = xyz           ; Condiciones periódicas all"""

    with open('ions.mdp', 'w') as f:
        f.write(content)
    print("✅ ions.mdp creado")

def create_em_mdp():
    """MDP para minimización de energía"""
    content = """; Archivo MDP para Minimización de Energía
; Sistema: Acetaminofen + PVA Polymer

integrator      = steep         ; Algoritmo steep descent
nsteps          = 50000         ; Máximo 50,000 pasos
emtol           = 10.0          ; Tolerancia fuerza (kJ/mol/nm)
emstep          = 0.01          ; Tamaño paso inicial (nm)

; Parámetros de salida
nstlog          = 500           ; Frecuencia log cada 500 pasos
nstenergy       = 500           ; Frecuencia energía cada 500 pasos

; Control temperatura - DESACTIVADO en minimización
tcoupl          = no            ; Sin termostato

; Control presión - DESACTIVADO en minimización
pcoupl          = no            ; Sin barostato

; Parámetros electrostáticos
coulombtype     = PME           ; Particle Mesh Ewald
rcoulomb        = 1.0           ; Radio corte coulomb (nm)
rvdw            = 1.0           ; Radio corte van der Waals (nm)
rlist           = 1.0           ; Radio lista vecinos (nm)

; Condiciones frontera
pbc             = xyz           ; Condiciones frontera periódicas

; Configuración PME
fourierspacing  = 0.16          ; Espaciado grilla PME (nm)
pme_order       = 4             ; Orden interpolación PME"""

    with open('em.mdp', 'w') as f:
        f.write(content)
    print("✅ em.mdp creado")

def create_nvt_mdp():
    """MDP para equilibración NVT (temperatura constante)"""
    content = """; Archivo MDP para Equilibración NVT
; Sistema: Acetaminofen + PVA Polymer
; Tiempo total: 100 ps

; Parámetros simulación
integrator      = md            ; Leap-frog integrator
nsteps          = 50000         ; 100 ps (50,000 * 0.002 ps)
dt              = 0.002         ; Tamaño paso 2 fs

; Parámetros salida
nstlog          = 1000          ; Log cada 2 ps
nstxout         = 1000          ; Coordenadas cada 2 ps
nstvout         = 1000          ; Velocidades cada 2 ps  
nstenergy       = 1000          ; Energía cada 2 ps
nstxout-compressed = 1000       ; Coordenadas comprimidas cada 2 ps

; Restricciones de enlace
constraints     = h-bonds       ; Restricción enlaces H
constraint_algorithm = lincs    ; Algoritmo LINCS
lincs_iter      = 1             ; Iteraciones LINCS
lincs_order     = 4             ; Orden LINCS

; Control temperatura - ACTIVADO
tcoupl          = V-rescale     ; Termostato V-rescale
tc-grps         = non-Water Water ; Grupos temperatura
tau_t           = 0.1    0.1    ; Constante tiempo (ps)
ref_t           = 300    300    ; Temperatura referencia (K)

; Control presión - DESACTIVADO en NVT
pcoupl          = no            ; Sin barostato

; Parámetros electrostáticos
coulombtype     = PME           ; Particle Mesh Ewald
rcoulomb        = 1.0           ; Radio corte coulomb (nm)
rvdw            = 1.0           ; Radio corte vdW (nm)
rlist           = 1.0           ; Radio lista vecinos (nm)

; Condiciones frontera
pbc             = xyz           ; Condiciones frontera periódicas

; Configuración PME
fourierspacing  = 0.16          ; Espaciado grilla PME (nm)

; Generación velocidades
gen_vel         = yes           ; Generar velocidades iniciales
gen_temp        = 300           ; Temperatura para velocidades (K)
gen_seed        = -1            ; Semilla aleatoria"""

    with open('nvt.mdp', 'w') as f:
        f.write(content)
    print("✅ nvt.mdp creado")

def create_npt_mdp():
    """MDP para equilibración NPT (presión y temperatura constantes)"""
    content = """; Archivo MDP para Equilibración NPT
; Sistema: Acetaminofen + PVA Polymer
; Tiempo total: 100 ps

; Parámetros simulación
integrator      = md            ; Leap-frog integrator
nsteps          = 50000         ; 100 ps (50,000 * 0.002 ps)
dt              = 0.002         ; Tamaño paso 2 fs

; Parámetros salida
nstlog          = 1000          ; Log cada 2 ps
nstxout         = 1000          ; Coordenadas cada 2 ps
nstvout         = 1000          ; Velocidades cada 2 ps
nstenergy       = 1000          ; Energía cada 2 ps
nstxout-compressed = 1000       ; Coordenadas comprimidas cada 2 ps

; Restricciones de enlace
constraints     = h-bonds       ; Restricción enlaces H
constraint_algorithm = lincs    ; Algoritmo LINCS
lincs_iter      = 1             ; Iteraciones LINCS
lincs_order     = 4             ; Orden LINCS

; Control temperatura - ACTIVADO
tcoupl          = V-rescale     ; Termostato V-rescale
tc-grps         = non-Water Water ; Grupos temperatura
tau_t           = 0.1    0.1    ; Constante tiempo (ps)
ref_t           = 300    300    ; Temperatura referencia (K)

; Control presión - ACTIVADO en NPT
pcoupl          = Parrinello-Rahman ; Barostato Parrinello-Rahman
pcoupltype      = isotropic     ; Acoplamiento isotrópico
tau_p           = 2.0           ; Constante tiempo presión (ps)
ref_p           = 1.0           ; Presión referencia (bar)
compressibility = 4.5e-5        ; Compresibilidad agua (bar^-1)

; Parámetros electrostáticos
coulombtype     = PME           ; Particle Mesh Ewald
rcoulomb        = 1.0           ; Radio corte coulomb (nm)
rvdw            = 1.0           ; Radio corte vdW (nm)
rlist           = 1.0           ; Radio lista vecinos (nm)

; Condiciones frontera
pbc             = xyz           ; Condiciones frontera periódicas

; Configuración PME
fourierspacing  = 0.16          ; Espaciado grilla PME (nm)

; Generación velocidades - NO en NPT
gen_vel         = no            ; Usar velocidades de NVT"""

    with open('npt.mdp', 'w') as f:
        f.write(content)
    print("✅ npt.mdp creado")

def create_md_mdp():
    """MDP para simulación de producción"""
    content = """; Archivo MDP para Simulación de Producción
; Sistema: Acetaminofen + PVA Polymer  
; Tiempo total: 10 ns

; Parámetros simulación
integrator      = md            ; Leap-frog integrator
nsteps          = 5000000       ; 10 ns (5,000,000 * 0.002 ps)
dt              = 0.002         ; Tamaño paso 2 fs

; Parámetros salida
nstlog          = 5000          ; Log cada 10 ps
nstxout         = 0             ; Sin coordenadas sin comprimir
nstvout         = 0             ; Sin velocidades
nstenergy       = 5000          ; Energía cada 10 ps
nstxout-compressed = 5000       ; Coordenadas comprimidas cada 10 ps
compressed-x-precision = 1000   ; Precisión coordenadas comprimidas

; Restricciones de enlace
constraints     = h-bonds       ; Restricción enlaces H
constraint_algorithm = lincs    ; Algoritmo LINCS
lincs_iter      = 1             ; Iteraciones LINCS
lincs_order     = 4             ; Orden LINCS

; Control temperatura
tcoupl          = V-rescale     ; Termostato V-rescale
tc-grps         = non-Water Water ; Grupos temperatura
tau_t           = 0.1    0.1    ; Constante tiempo (ps)
ref_t           = 300    300    ; Temperatura referencia (K)

; Control presión
pcoupl          = Parrinello-Rahman ; Barostato Parrinello-Rahman
pcoupltype      = isotropic     ; Acoplamiento isotrópico
tau_p           = 2.0           ; Constante tiempo presión (ps)
ref_p           = 1.0           ; Presión referencia (bar)
compressibility = 4.5e-5        ; Compresibilidad agua (bar^-1)

; Parámetros electrostáticos
coulombtype     = PME           ; Particle Mesh Ewald
rcoulomb        = 1.0           ; Radio corte coulomb (nm)
rvdw            = 1.0           ; Radio corte vdW (nm)
rlist           = 1.0           ; Radio lista vecinos (nm)

; Condiciones frontera
pbc             = xyz           ; Condiciones frontera periódicas

; Configuración PME
fourierspacing  = 0.16          ; Espaciado grilla PME (nm)

; Generación velocidades - NO en producción
gen_vel         = no            ; Usar velocidades de NPT

; Configuración avanzada para eficiencia
nstlist         = 10            ; Frecuencia actualización lista vecinos
ns_type         = grid          ; Búsqueda vecinos tipo grilla
cutoff-scheme   = Verlet        ; Esquema corte Verlet"""

    with open('md.mdp', 'w') as f:
        f.write(content)
    print("✅ md.mdp creado")

def verify_mdp_files():
    """Verificar archivos .mdp creados"""
    mdp_files = ['temporal.mdp', 'ions.mdp', 'em.mdp', 'nvt.mdp', 'npt.mdp', 'md.mdp']
    
    print("\n📋 VERIFICACIÓN DE ARCHIVOS .MDP:")
    print("=" * 35)
    
    for mdp_file in mdp_files:
        try:
            with open(mdp_file, 'r') as f:
                lines = f.readlines()
            
            # Buscar parámetros clave
            integrator = next((line for line in lines if 'integrator' in line and '=' in line), 'No encontrado')
            nsteps = next((line for line in lines if 'nsteps' in line and '=' in line), 'No encontrado')
            
            print(f"📄 {mdp_file}:")
            print(f"  Líneas: {len(lines)}")
            print(f"  Integrator: {integrator.split('=')[1].split(';')[0].strip() if '=' in integrator else 'No encontrado'}")
            print(f"  Nsteps: {nsteps.split('=')[1].split(';')[0].strip() if '=' in nsteps else 'No encontrado'}")
            print("  ✅ Archivo creado correctamente")
            
        except FileNotFoundError:
            print(f"❌ {mdp_file} NO CREADO")
        except Exception as e:
            print(f"⚠️ {mdp_file} - Error: {e}")
        
        print()

if __name__ == "__main__":
    print("🔧 GENERANDO TODOS LOS ARCHIVOS .MDP")
    print("=" * 40)
    print()
    
    # Crear todos los archivos .mdp
    create_temporal_mdp()
    create_ions_mdp()
    create_em_mdp()
    create_nvt_mdp()
    create_npt_mdp()
    create_md_mdp()
    
    # Verificar archivos creados
    verify_mdp_files()
    
    print("🎯 ARCHIVOS .MDP GENERADOS EXITOSAMENTE")
    print("=" * 45)
    print("📋 FLUJO COMPLETO DE SIMULACIÓN:")
    print("  1. temporal.mdp  → Verificación topología")
    print("  2. ions.mdp      → Añadir iones")
    print("  3. em.mdp        → Minimización energía") 
    print("  4. nvt.mdp       → Equilibración temperatura")
    print("  5. npt.mdp       → Equilibración presión")
    print("  6. md.mdp        → Simulación producción (10 ns)")
    print()
    print("✅ LISTO PARA CONTINUAR CON GENERACIÓN DE CAJA Y SOLVATACIÓN")
