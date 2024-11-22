import fitz  # PyMuPDF

with fitz.open('chall.pdf') as doc:
    for page in doc:
        print(page.get_text("text"))