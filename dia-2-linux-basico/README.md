# ▶ Día 2: Linux Básico y ACPYPE

**Orden de ejecución:**
1. **[paso-1-estructura-directorios.sh](paso-1-estructura-directorios.sh)** → Crear estructura de carpetas del proyecto
2. **[paso-2-verificar-sistema.sh](paso-2-verificar-sistema.sh)** → Verificar arquitectura y compatibilidad
3. **[paso-3-instalar-miniconda.sh](paso-3-instalar-miniconda.sh)** → Instalar Miniconda para gestión de paquetes
4. **[paso-4-instalar-ambertools.sh](paso-4-instalar-ambertools.sh)** → Instalar AmberTools (incluye ACPYPE)
5. **[paso-5-procesar-moleculas.sh](paso-5-procesar-moleculas.sh)** → Procesar moléculas .mol2 con ACPYPE
6. **[paso-6-organizar-resultados.sh](paso-6-organizar-resultados.sh)** → Organizar archivos para GROMACS

## ■ Dónde ejecutar cada paso:
- **Todos los pasos:** Ubuntu (WSL2)

## ⚠ Importante:
- Tener archivos .mol2 listos antes del paso 5
- El paso 3 requiere reiniciar terminal (source ~/.bashrc)
- El paso 4 puede tomar varios minutos
- No interrumpir el procesamiento ACPYPE en paso 5
