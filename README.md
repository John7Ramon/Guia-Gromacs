<div align="center">

# Repositorio Guía de Comandos GROMACS


**Banco completo de comandos para acompañar la guía práctica de GROMACS en Windows**

![GROMACS](https://img.shields.io/badge/GROMACS-2023-2E86AB?style=flat-square)
![WSL](https://img.shields.io/badge/WSL-Ubuntu-E95420?style=flat-square)
![Python](https://img.shields.io/badge/Python-3.8+-3776AB?style=flat-square)
![Licencia](https://img.shields.io/badge/Licencia-MIT-green?style=flat-square)

*Repositorio estructurado diseñado para estudiantes y profesionales que aprenden simulaciones de dinámica molecular con GROMACS*

</div>

---

## Descripción General

Este repositorio proporciona una referencia completa de comandos diseñada para acompañar la guía práctica de GROMACS. Ofrece comandos paso a paso organizados por módulos de aprendizaje, permitiendo a los usuarios progresar desde la instalación básica hasta simulaciones avanzadas de dinámica molecular.

### Flujo de Trabajo

1. **Integración de Referencias**: Consulte su guía PDF para referencias de comandos
2. **Navegación por Módulos**: Acceda a la carpeta del día correspondiente
3. **Ejecución de Comandos**: Copie y ejecute comandos de archivos organizados
4. **Verificación**: Confirme que los resultados coincidan con las salidas esperadas
5. **Seguimiento del Progreso**: Siga la ruta de aprendizaje estructurada

---

## Estructura del Repositorio

### Día 1: Instalación de WSL y GROMACS
**Estado**: Completo | **Plataforma**: Windows + WSL | **Duración**: ~2 horas

| Componente | Descripción | Entorno | Tiempo |
|------------|-------------|---------|---------|
| Preparación del Sistema | Actualizaciones de BIOS, Windows, Visual Studio | Windows | 30 min |
| Comandos PowerShell | Comandos DISM para habilitar WSL | PowerShell (Admin) | 10 min |
| Instalación WSL | Configuración WSL post-reinicio | PowerShell (Admin) | 20 min |
| Instalación GROMACS | Siguiendo tutorial específico de YouTube | Ubuntu WSL | 45 min |
| Verificación de Instalación | Verificar funcionalidad de GROMACS | Ubuntu WSL | 15 min |

**Resultado Esperado**: Entorno WSL funcional con instalación verificada de GROMACS

### Día 2: Fundamentos de Linux y ACPYPE
**Estado**: Completo | **Plataforma**: Ubuntu WSL | **Duración**: ~1.5 horas

| Componente | Descripción | Herramientas | Tiempo |
|------------|-------------|--------------|--------|
| Estructura de Directorios | Organización de carpetas del proyecto | `mkdir`, `ls` | 10 min |
| Verificación del Sistema | Verificación de arquitectura y compatibilidad | `uname`, `lscpu` | 5 min |
| Instalación Miniconda | Configuración del sistema de gestión de paquetes | `wget`, `bash` | 15 min |
| Instalación AmberTools | ACPYPE y herramientas moleculares | `conda` | 30 min |
| Procesamiento Molecular | Procesamiento de archivos .mol2 con ACPYPE | `acpype` | 20 min |
| Organización de Archivos | Preparación de archivos listos para GROMACS | `mv`, `cp` | 10 min |

**Resultado Esperado**: Entorno Python científico con capacidades de procesamiento molecular

### Día 3: Preparación de Sistemas Moleculares
**Estado**: Completo | **Plataforma**: Ubuntu WSL | **Duración**: ~3 horas

| Componente | Descripción | Herramientas | Complejidad |
|------------|-------------|--------------|-------------|
| Validación de Archivos | Verificación de consistencia molecular | `gmx check`, scripts personalizados | Intermedio |
| Diagnóstico del Sistema | Análisis de estructura y compatibilidad | `gmx editconf`, análisis manual | Intermedio |
| Consolidación de Atomtypes | Resolución de atomtypes duplicados | `awk`, `sed`, scripts bash | Avanzado |
| Adición de Agua | Integración de moléculas de agua | `gmx insert-molecules` | Intermedio |
| Generación del Sistema | Creación completa del sistema molecular | Python, MDAnalysis | Avanzado |
| Verificación de Topología | Validación del sistema generado | `gmx grompp` | Intermedio |
| Generación de Parámetros | Creación de archivos de parámetros MD | Python, plantillas | Avanzado |
| Caja de Simulación | Definición de geometría de caja | `gmx editconf` | Intermedio |
| Solvatación del Sistema | Adición de entorno acuoso | `gmx solvate` | Intermedio |
| Reparación de Topología | Correcciones de topología post-solvatación | Edición manual, scripts | Avanzado |
| Neutralización con Iones | Neutralización de carga del sistema | `gmx genion` | Intermedio |
| Minimización de Energía | Optimización del sistema final | `gmx mdrun` | Intermedio |

**Resultado Esperado**: Sistema molecular completo listo para simulaciones MD

### Hoja de Ruta de Desarrollo

| Módulo | Contenido | Estado | Lanzamiento Estimado |
|--------|-----------|---------|---------------------|
| **Día 4** | Ejecución de Simulación MD | En Desarrollo | Q4 2025 |
| **Día 5** | Análisis de Trayectorias | Planificado | Q1 2026 |
| **Día 6** | Visualización Avanzada | Planificado | Q1 2026 |

---

## Requisitos del Sistema

### Especificaciones de Hardware
- **Memoria**: 8 GB RAM mínimo (16 GB recomendado)
- **Almacenamiento**: 50 GB de espacio disponible
- **Procesador**: CPU de cuatro núcleos (8+ núcleos recomendado)
- **Gráficos**: GPU opcional (aceleración significativa para simulaciones grandes)

### Dependencias de Software
- **Sistema Operativo**: Windows 10/11 (build 1903+)
- **WSL**: Versión 2 habilitada
- **Shell**: PowerShell con privilegios administrativos
- **Red**: Conexión a internet para descargas de paquetes

### Verificación del Entorno

```bash
# Comandos de verificación del sistema
echo "Sistema: $(uname -a)"
echo "GROMACS: $(gmx --version 2>/dev/null | head -n1 || echo 'No instalado')"
echo "Python: $(python3 --version 2>/dev/null || echo 'No instalado')"
```

### Requisitos de Plataforma por Módulo

| Herramienta | Día 1 (Pasos 1-3) | Día 1 (Pasos 4-5) | Día 2+ |
|-------------|--------------------|--------------------|---------|
| **PowerShell (Admin)** | Requerido | No usado | No usado |
| **Ubuntu WSL** | No usado | Requerido | Requerido |
| **GROMACS** | No disponible | Se instala | Requerido |
| **Python/Conda** | No disponible | No usado | Se instala |

---

## Solución de Problemas

### Problemas de Instalación WSL

```powershell
# Verificar soporte de virtualización
systeminfo | find "Hyper-V"

# Habilitar características requeridas de Windows
dism /online /enable-feature /featurename:Microsoft-Windows-Subsystem-Linux
dism /online /enable-feature /featurename:VirtualMachinePlatform
```

### Comando GROMACS No Encontrado

```bash

# Recargar variables de entorno
source ~/.bashrc
```

### Problemas del Entorno Conda

```bash
# Verificar instalación de conda
ls ~/miniconda3/bin/conda

```

---

## Contribuciones

Las contribuciones son bienvenidas a través del flujo de trabajo estándar de GitHub:

1. Hacer fork del repositorio
2. Crear una rama de funcionalidad (`git checkout -b feature/mejora`)
3. Confirmar cambios (`git commit -am 'Agregar nueva funcionalidad'`)
4. Push a la rama (`git push origin feature/mejora`)
5. Crear un Pull Request

### Directrices para Contribuciones
- Seguir el estilo de código y formato de documentación existente
- Incluir mensajes de commit claros
- Probar comandos en instalación WSL limpia
- Actualizar documentación para nuevas funcionalidades

---

## Información de Contacto

**John Byron Ramón Pérez**  

- **Correo**: jbramon2@utpl.edu.ec
- **Institución**: Universidad Técnica Particular de Loja (UTPL)

Para consultas técnicas, reportes de errores o solicitudes de funcionalidades, utilice el sistema de seguimiento de issues del repositorio o contacte por correo electrónico.

---


**Desarrollado para la comunidad de simulación molecular**

*Si encuentra útil este repositorio, considere darle una estrella para ayudar a otros a descubrir estos recursos*

</div>
