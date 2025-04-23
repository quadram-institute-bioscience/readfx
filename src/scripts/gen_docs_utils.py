#!/usr/bin/env python3

import re
import os

def extract_methods_with_docs(file_path):
    """
    Extract method declarations and docstrings from Nim files
    
    Args:
        file_path: Path to the Nim file
        
    Returns:
        List of dictionaries with method name, declaration, and docstring
    """
    with open(file_path, 'r') as f:
        content = f.read()
    
    # Find all proc declarations and their docstrings
    pattern = r'proc\s+(\w+)(\*?)\s*\((.*?)\)(?:\s*:\s*(\w+))?\s*=\s*\n((?:\s*##.*\n)+)'
    methods = []
    
    # Find all regular proc declarations with docstrings
    matches = re.finditer(pattern, content)
    for match in matches:
        name = match.group(1)
        is_exported = bool(match.group(2))
        params = match.group(3)
        return_type = match.group(4) or ""
        docstring = match.group(5)
        
        # Clean up the docstring
        doc_lines = [line.strip().replace('## ', '').replace('##', '')
                    for line in docstring.split('\n') if '##' in line]
        doc_text = '\n'.join(doc_lines)
        
        # Format the declaration
        if return_type:
            declaration = f"proc {name}{match.group(2)}({params}): {return_type} ="
        else:
            declaration = f"proc {name}{match.group(2)}({params}) ="
            
        methods.append({
            'name': name,
            'declaration': declaration,
            'docstring': doc_text,
            'is_exported': is_exported
        })
    
    # Handle other proc declarations with different patterns
    alt_pattern = r'proc\s+(\w+)(\*?)\s*\((.*?)\)(?:\s*:\s*(\w+))?\s*=\s*[^\n]+\n((?:\s*##.*\n)+)'
    matches = re.finditer(alt_pattern, content)
    for match in matches:
        # Skip if this method was already found
        name = match.group(1)
        if any(m['name'] == name for m in methods):
            continue
            
        is_exported = bool(match.group(2))
        params = match.group(3)
        return_type = match.group(4) or ""
        docstring = match.group(5)
        
        # Clean up the docstring
        doc_lines = [line.strip().replace('## ', '').replace('##', '')
                    for line in docstring.split('\n') if '##' in line]
        doc_text = '\n'.join(doc_lines)
        
        # Format the declaration
        if return_type:
            declaration = f"proc {name}{match.group(2)}({params}): {return_type} ="
        else:
            declaration = f"proc {name}{match.group(2)}({params}) ="
            
        methods.append({
            'name': name,
            'declaration': declaration,
            'docstring': doc_text,
            'is_exported': is_exported
        })
        
    return methods

def generate_markdown(methods, output_file):
    """
    Generate markdown documentation from extracted methods
    
    Args:
        methods: List of method dictionaries
        output_file: Path to output markdown file
    """
    with open(output_file, 'w') as f:
        f.write("# ReadFX Utility Functions\n\n")
        f.write("This file documents the utility functions available in the ReadFX library.\n\n")
        
        for method in methods:
            name = method['name']
            name_with_export = f"{name}*" if method['is_exported'] else name
            
            f.write(f"## {name_with_export}\n\n")
            f.write("```nim\n")
            f.write(f"{method['declaration']}\n")
            f.write("```\n\n")
            f.write(f"{method['docstring']}\n\n")

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, '..', '..'))
    
    sequtils_path = os.path.join(project_root, 'readfx', 'sequtils.nim')
    output_path = os.path.join(project_root, 'docs', 'UTILS.md')
    
    methods = extract_methods_with_docs(sequtils_path)
    generate_markdown(methods, output_path)
    
    print(f"Documentation generated: {output_path}")
    print(f"Extracted {len(methods)} methods from sequtils.nim")

if __name__ == "__main__":
    main()