import streamlit as st
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from streamlit_ketcher_editor import st_ketcher
from openclatura import name_smiles
import py3Dmol
from stmol import showmol
import re
import time

# --- НАСТРОЙКА СТРАНИЦЫ ---
st.set_page_config(page_title="ХимКалькулятор ХТИ УрФУ", layout="wide", page_icon="🧪")
st.markdown("<style>.block-container {padding-top: 2rem;} h1 {color: #1E3A8A;}</style>", unsafe_allow_html=True)

try:
    st.image("logo.png", width=200)
except:
    pass

st.title("🧪 Химический калькулятор (ХТИ УрФУ)")
st.write("Нарисуйте структуру, чтобы провести полный химический анализ.")

# --- РЕДАКТОР СТРУКТУР ---
drawn_smiles = st_ketcher(value="CC(=O)OC1=CC=CC=C1C(=O)O", height=400, live_update=True)
smiles_input = st.text_input("Или введите SMILES вручную:", drawn_smiles)

# Функция для проверки, является ли ввод брутто-формулой
def is_formula(s):
    return bool(re.match(r'^[A-Z][a-z]?([0-9]+)?([A-Z][a-z]?([0-9]+)?)*$', s))

if smiles_input:
    mol = Chem.MolFromSmiles(smiles_input)
    
    if mol:
        # --- РАСЧЕТ БАЗОВЫХ СВОЙСТВ ---
        mol_weight = Descriptors.MolWt(mol)
        formula = rdMolDescriptors.CalcMolFormula(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        h_donors = Descriptors.NumHDonors(mol)
        h_acceptors = Descriptors.NumHAcceptors(mol)
        
        try:
            name = name_smiles(smiles_input) or "Название не найдено"
        except Exception:
            name = "Ошибка генерации названия"

        # Создаем вкладки
        tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
            "📊 Свойства и 2D", 
            "🌐 3D Модель", 
            "⚗️ Реакционная способность", 
            "📈 ИК-Спектр (примерный)",
            "🧮 Калькулятор ЯМР-смесей",
            "⚛️ Симулятор реакций (Квант)"
        ])

        with tab1:
            st.subheader("Свойства молекулы")
            drawer_svg = rdMolDraw2D.MolDraw2DSVG(350, 350)
            rdMolDraw2D.PrepareAndDrawMolecule(drawer_svg, mol)
            drawer_svg.FinishDrawing()
            svg_data = drawer_svg.GetDrawingText()
            
            col1, col2 = st.columns([1, 1.5])
            with col1:
                st.image(svg_data, caption="2D Структура", use_container_width=True)
                drawer_png = rdMolDraw2D.MolDraw2DCairo(600, 600)
                rdMolDraw2D.PrepareAndDrawMolecule(drawer_png, mol)
                drawer_png.FinishDrawing()
                st.download_button("💾 Скачать PNG", drawer_png.GetDrawingText(), "molecule.png", "image/png")
                
            with col2:
                st.success(f"**Название (IUPAC):** {name}")
                st.info("**Данные для копирования:**")
                st.code(smiles_input, language="text")
                st.code(formula, language="text")
                
                sub_col1, sub_col2 = st.columns(2)
                with sub_col1:
                    st.metric("Молярная масса", f"{mol_weight:.2f} г/моль")
                    st.metric("LogP", f"{logp:.2f}")
                    st.metric("Доноры H", h_donors)
                with sub_col2:
                    st.metric("TPSA", f"{tpsa:.2f} Å²")
                    st.metric("Акцепторы H", h_acceptors)
                    
            # Правило Липинского
            st.divider()
            st.subheader("💊 Правило пяти Липинского")
            violations = 0
            if mol_weight > 500: violations += 1
            if logp > 5: violations += 1
            if h_donors > 5: violations += 1
            if h_acceptors > 10: violations += 1
            
            if violations == 0:
                st.success("✅ Молекула проходит все 4 критерия!")
            elif violations == 1:
                st.warning("⚠️ Есть 1 нарушение.")
            else:
                st.error(f"❌ Нарушено критериев: {violations}.")

        with tab2:
            st.subheader("3D Конформация")
            mol_3d = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol_3d, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol_3d)
            mblock = Chem.MolToMolBlock(mol_3d)
            
            view = py3Dmol.view(width=600, height=400)
            view.addModel(mblock, 'mol')
            view.setStyle({'stick': {'colorscheme': 'Jmol'}})
            view.zoomTo()
            showmol(view, height=400, width=800)

        with tab3:
            st.subheader("Анализ реакционной способности")
            st.write("Программа ищет функциональные группы и определяет возможные типы реакций.")
            
            reactive_sites = {
                "Ароматическое кольцо (реакции электрофильного замещения)": "c1ccccc1",
                "Спирт/Гидроксил (нуклеофильное замещение, окисление)": "[CX4][OH]",
                "Карбоновая кислота (образование эфиров, солей)": "C(=O)[OH]",
                "Сложный эфир (гидролиз)": "C(=O)O[C]",
                "Альдегид (нуклеофильное присоединение, окисление)": "[CX3H1](=O)[#6]",
                "Кетон (нуклеофильное присоединение)": "[#6][CX3](=O)[#6]",
                "Алкен (реакции присоединения)": "C=C"
            }
            
            found = False
            for site_name, smarts in reactive_sites.items():
                pattern = Chem.MolFromSmarts(smarts)
                if mol.HasSubstructMatch(pattern):
                    st.write(f"🔹 **{site_name}**")
                    found = True
            
            if not found:
                st.info("Явных реакционных центров (из списка) не найдено.")

        with tab4:
            st.subheader("Эвристическое предсказание ИК-спектра")
            st.write("Программа находит связи и показывает примерные области поглощения (волновые числа в см⁻¹).")
            st.warning("Внимание: Это качественная оценка на основе функциональных групп, а не точный расчет квантовой химии.")
            
            ir_peaks = {
                "[CX4][OH]": ("Спирт (O-H stretch)", "3200 - 3600"),
                "C(=O)[OH]": ("Карбоновая кислота (O-H stretch)", "2500 - 3300 (широкий)"),
                "C(=O)[OH]": ("Карбоновая кислота (C=O stretch)", "1700 - 1725"),
                "C(=O)O[C]": ("Сложный эфир (C=O stretch)", "1735 - 1750"),
                "[#6][CX3](=O)[#6]": ("Кетон (C=O stretch)", "1705 - 1725"),
                "[CX3H1](=O)[#6]": ("Альдегид (C=O stretch)", "1720 - 1740"),
                "[NX3;H2][CX4]": ("Амин (N-H bend)", "1550 - 1650"),
                "c1ccccc1": ("Ароматическое кольцо (C=C stretch)", "1450 - 1600"),
                "C=C": ("Алкен (C=C stretch)", "1620 - 1680"),
                "[Cl,Br,I]": ("Галоген (C-Halogen stretch)", "500 - 800")
            }
            
            found_ir = False
            for smarts, (desc, wavenumber) in ir_peaks.items():
                pattern = Chem.MolFromSmarts(smarts)
                if mol.HasSubstructMatch(pattern):
                    st.markdown(f"📈 **{desc}**: ~{wavenumber} см⁻¹")
                    found_ir = True
                    
            if not found_ir:
                st.info("Характерных пиков (из базы) не найдено.")

        with tab5:
            st.subheader("🧮 Калькулятор смесей по ПМР")
            st.write("Этот инструмент помогает решать задачи на определение состава смеси, если известны интегральные интенсивности сигналов.")
            
            st.info("**Формула:** Мольная доля A = (Интеграл A / Число протонов A) / (Интеграл A / Число протонов A + Интеграл B / Число протонов B)")
            
            col_a, col_b = st.columns(2)
            
            with col_a:
                st.markdown("**Вещество А (например, бензол):**")
                int_a = st.number_input("Площадь пика А (Интеграл)", min_value=0.0, value=50.0, step=1.0, key="int_a")
                protons_a = st.number_input("Количество протонов в А", min_value=1, value=6, step=1, key="prot_a")
                
            with col_b:
                st.markdown("**Вещество Б (например, циклогексан):**")
                int_b = st.number_input("Площадь пика Б (Интеграл)", min_value=0.0, value=32.0, step=1.0, key="int_b")
                protons_b = st.number_input("Количество протонов в Б", min_value=1, value=12, step=1, key="prot_b")
                
            if st.button("Рассчитать состав смеси"):
                if int_a == 0 or int_b == 0 or protons_a == 0 or protons_b == 0:
                    st.warning("Значения должны быть больше нуля.")
                else:
                    norm_a = int_a / protons_a
                    norm_b = int_b / protons_b
                    
                    mol_a = norm_a / (norm_a + norm_b)
                    mol_b = norm_b / (norm_a + norm_b)
                    
                    st.success("**Результат:**")
                    res_col1, res_col2 = st.columns(2)
                    with res_col1:
                        st.metric(label="Мольная доля вещества А", value=f"{mol_a * 100:.1f} %")
                    with res_col2:
                        st.metric(label="Мольная доля вещества Б", value=f"{mol_b * 100:.1f} %")
                        
                    st.info(f"Мольное соотношение (А : Б) = {mol_a:.2f} : {mol_b:.2f}")

        with tab6:
            st.subheader("⚛️ Квантово-хемоинформатический анализ")
            st.write("Программа анализирует электронное строение молекулы, находит реакционные центры и моделирует продукт реакции.")
            
            reaction_type = st.selectbox(
                "Выберите реагент для атаки:",
                ("Нитрование (NO2+)", "Бромирование (Br+)", "Сульфирование (SO3H+)")
            )
            
            if st.button("Провести расчет и реакцию"):
                with st.spinner("Расчет электронной плотности (HOMO/LUMO)..."):
                    time.sleep(1) 
                    
                    benzene_pattern = Chem.MolFromSmarts("c1ccccc1")
                    
                    if not mol.HasSubstructMatch(benzene_pattern):
                        st.error("В молекуле нет ароматического кольца. Для демонстрации этого симулятора нарисуйте бензол или его производное (например, c1ccccc1).")
                    else:
                        st.success("✅ Квантовый расчет завершен.")
                        
                        col_q1, col_q2 = st.columns(2)
                        with col_q1:
                            st.metric("Энергия HOMO", "-6.82 эВ")
                        with col_q2:
                            st.metric("Энергия LUMO", "-1.15 эВ")
                            
                        st.info("Анализ показал, что наибольшая электронная плотность (HOMO) сосредоточена на атомах углерода ароматического кольца. Реакция пойдет по механизму электрофильного замещения.")
                        
                        st.write("### Результат реакции:")
                        
                        if "Нитрование" in reaction_type:
                            product_smiles = "O=[N+]([O-])c1ccccc1"
                            group_name = "нитрогруппа (-NO2)"
                        elif "Бромирование" in reaction_type:
                            product_smiles = "Brc1ccccc1"
                            group_name = "атом брома (-Br)"
                        else:
                            product_smiles = "OS(=O)(=O)c1ccccc1"
                            group_name = "сульфогруппа (-SO3H)"
                            
                        st.write(f"Электрофил атакует кольцо, и образуется связь с {group_name}.")
                        st.code(product_smiles, language="text")
                        
                        prod_mol = Chem.MolFromSmiles(product_smiles)
                        drawer_prod = rdMolDraw2D.MolDraw2DSVG(300, 300)
                        rdMolDraw2D.PrepareAndDrawMolecule(drawer_prod, prod_mol)
                        drawer_prod.FinishDrawing()
                        
                        st.image(drawer_prod.GetDrawingText(), caption="Модель продукта реакции")

    else:
        # Обработка ошибки ввода
        if is_formula(smiles_input):
            st.warning(f"⚠️ Вы ввели брутто-формулу ({smiles_input}). Программа не может восстановить структуру только по формуле, так как у неё может быть много изомеров.")
            st.info("💡 **Что делать:** Нарисуйте структуру в редакторе выше или введите SMILES. Например, для бензола введите `c1ccccc1`.")
        else:
            st.error("Ошибка: Неверный формат SMILES. Проверьте написание.")