# ▶ Día 3: Preparación de Sistemas Moleculares

**Orden de ejecución:**

1. **[01-validacion-archivos.txt](01-validacion-archivos.txt)** → Validación completa de consistencia molecular
2. **[02-diagnostico-completo.txt](02-diagnostico-completo.txt)** → Análisis profundo de estructura y compatibilidad
3. **[03-consolidacion-atomtypes.txt](03-consolidacion-atomtypes.txt)** → Consolidar atomtypes duplicados en archivo único
4. **[04-agregar-agua-archivos-limpios.txt](04-agregar-agua-archivos-limpios.txt)** → Agregar moléculas de agua al sistema
5. **[05-crear-sistema.py](05-crear-sistema.py)** → Script para generar sistema molecular completo
6. **[06-verificación-topologia.txt](06-verificación-topologia.txt)** → Verificar topología del sistema generado
7. **[07-generar-mdp.py](07-generar-mdp.py)** → Script para generar archivos de parámetros MD
8. **[08-generar-caja.txt](08-generar-caja.txt)** → Crear caja de simulación para el sistema
9. **[09-solvatar-sistema.txt](09-solvatar-sistema.txt)** → Solvatar sistema con agua
10. **[10-reparar-topologia-agua.txt](10-reparar-topologia-agua.txt)** → Reparar topología después de solvatación  
11. **[11-neutralizar-iones.txt](11-neutralizar-iones.txt)** → Neutralizar sistema con iones
12. **[12-minimizacion-energia.txt](12-minimizacion-energia.txt)** → Minimizar energía del sistema final

## ■ Dónde ejecutar cada paso:
- **Todos los pasos:** Ubuntu (WSL2) con GROMACS y AmberTools instalados

## ⚠ Importante:
- **Paso 1:** Requiere archivos .gro/.itp existentes para validar
- **Paso 2:** Analiza estructura interna para detectar conflictos de atomtypes (puede tomar tiempo)
- **Paso 5 y 7:** Hacer ejecutables los scripts Python: `python3 nombre_script.py`
- **Pasos 8-12:** Proceso secuencial de preparación del sistema - seguir orden exacto
- **Paso 12:** La minimización puede tomar varios minutos dependiendo del tamaño del sistema
