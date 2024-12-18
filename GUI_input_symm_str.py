################################################################################################################
# python GUI interface to create the input file to run the fortran symmetric structure generation symm_str code#
################################################################################################################

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
import re
import os
#import program_main
import subprocess
import sys
import cisvb
import numpy as np


class Ctrl_Input:
    def __init__(self, root):
# initialise disctionary to store control data
        self.ctrl_inputs={}
        self.input_text=''
        self.fname=''
        self.other_inp=True
 
        # Main Frame
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Input Fields
        self.entries = {}
        self.keywords = ["str", "nao", "nae", "nmul", "frgtyp", "chinst"]
        self.create_ctrl_pans()



    def create_ctrl_pans(self):
        ttk.Label(self.frame, text="Enter Ctrl Keywords").grid(row=1, column=0, columnspan=2, pady=5)
        ttk.Label(self.frame, text="file_name").grid(row=0, column=0, sticky=tk.W, padx=5, pady=5)
        self.file_name = ttk.Entry(self.frame, width=20)
        self.file_name.grid(row=0, column=1, padx=5, pady=5)
#        self.entries["file_name"]=entry
        for idx, key in enumerate(self.keywords):
            ttk.Label(self.frame, text=key).grid(row=idx + 2, column=0, sticky=tk.W, padx=5, pady=5)
            entry = ttk.Entry(self.frame, width=20)
            entry.grid(row=idx + 2, column=1, padx=5, pady=5)
            self.entries[key] = entry



#    def go_to_output(self):
#        self.frame.destroy()  # Destroy the output frame
#        InputCreator(self.root)  # Switch back to the input frame
#

#    def reset_fields(self):
#        for entry in self.entries.values():
#            entry.delete(0, tk.END)
#        if self.fragment_entries:
#            for pane in self.fragment_entries:
#                for field in pane:
#                    field.delete(0, tk.END)
#        if self.orbital_entries:
#            for pane in self.orbital_entries:
#                for field in pane:
#                    field.delete(0, tk.END)

###################################################################################
########## orbital section starts here :
###################################################################################

class Orb_Input:
    def __init__(self,root):
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        self.orbital_frame = None
        self.assoatm_entries = []
        self.assotyp_entries = []
        self.description_orb = ""


    def create_orbital_section(self):
        # Create a new frame for orbital inputs if it doesn't exist
        if self.orbital_frame is None:
            self.orbital_frame = tk.Toplevel(self.root, padx=10, pady=10)
            self.orbital_frame.title("orbital inputs")
            self.orbital_frame.geometry("600x600")

        # Clear existing fragment panes
        for widget in self.orbital_frame.winfo_children():
            widget.destroy()
        self.orbital_entries = []
        self.assoatm_entries = []
        self.assotyp_entries = []
        self.analyze_orbital()


    def analyze_orbital(self):
        self.other_inp=False
        ctl_key = self.generate_input()  # calling generate_input to get ctrl values
        num_orbital = int(ctl_key["nao"])  # read the number of active orbitals 
#        print('nao', num_orbital)
        self.other_inp=True

        # Clear existing panes
        for widget in self.orbital_frame.winfo_children():
            widget.destroy()

        # Create input panes for each active orbitals
        ttk.Label(self.orbital_frame, text=f"Number of active orbitals: {num_orbital}").grid(
            row=0, column=0, columnspan=2, pady=5
        )

#       creating the pane for getting orbital numbers

        ttk.Label(self.orbital_frame, text=f"associated atom num").grid(row=2, column=3, padx=5, pady=5, sticky=tk.W)
#        ttk.Label(self.orbital_frame, text=f" Orbital type").grid(row=1, column=4, padx=5, pady=1, sticky=tk.W)
        ttk.Label(self.orbital_frame, text=f"Sig-Pi/fragment").grid(row=2, column=4, padx=5, pady=5, sticky=tk.W)
        for i in range(num_orbital):
            ttk.Label(self.orbital_frame, text=f"active orbital {i+1} ").grid(row=i+3, column=1, padx=5, pady=5, sticky=tk.W)
#       creating the pane for getting associated atom number
            assoatm_entry = ttk.Entry(self.orbital_frame, width=10)
            assoatm_entry.grid(row=i+3, column=3, padx=5, pady=5)

            self.assoatm_entries.append(assoatm_entry)
#       creating the pane for getting associated atom number
            assotyp_entry = ttk.Entry(self.orbital_frame, width=10)
            assotyp_entry.grid(row=i+3, column=4, padx=5, pady=5)

            self.assotyp_entries.append(assotyp_entry)
        print(self.assotyp_entries, self.assoatm_entries)


    def calculate_orbital(self, description_orb):
        """
        Parse the description and calculate the total number of fragments.
        """
        total_orbital = 0
        groups = description_orb.split()
        for group in groups:
            if '*' in group:
                # Handle compact notation like "3*4"
                match = re.match(r"(\d+)\*(\d+)", group)
                if match:
                    total_orbital += int(match.group(2))
            else:
                # Handle individual groups
                total_orbital += 1
        return total_orbital


################################################################################
#### Spatial keyword section starts here : first ##
################################################################################

class Spl_Keywd:
    def __init__(self, root):
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.keywd_frame = None
        self.keywd_entries = []
        self.description_key = ""

    def create_keywd_section(self):
        # Create a new frame for fragment inputs if it doesn't exist
        if self.keywd_frame is None:
            self.keywd_frame = tk.Toplevel(self.root, padx=10, pady=10)
            self.keywd_frame.title("Spatial Keywords")
            self.keywd_frame.geometry("400x300")

        # Input for description of fragments
        ttk.Label(self.keywd_frame, text="Number of Spatial Inputs").grid(row=0, column=0, columnspan=2, pady=5)
        desc_entry = ttk.Entry(self.keywd_frame, width=40)
        desc_entry.grid(row=1, column=0, columnspan=2, pady=5)
        ttk.Button(
            self.keywd_frame, text="Create Inputs", command=lambda: self.analyze_keywds(desc_entry.get())
        ).grid(row=2, column=0, columnspan=2, pady=10)

    def analyze_keywds(self, description):
        self.description_key = int(description)
        try:
            # Calculate the total number of fragments
            num_keywds = self.description_key
            if num_keywds <= 0:
                raise ValueError("Number of keywords must be positive.")
        except ValueError:
            messagebox.showerror("Invalid Input", "Please enter a valid number of keywords.")
            return

        # Clear existing key word panes
        for widget in self.keywd_frame.winfo_children():
            widget.destroy()

        # Create input panes for each fragment
        ttk.Label(self.keywd_frame, text=f"Number of keywords: {num_keywds}").grid(
            row=0, column=0, columnspan=2, pady=5
        )
        for i in range(num_keywds):
            row = i + 1
            ttk.Label(self.keywd_frame, text=f"keywords {i + 1} Type").grid(row=row, column=0, padx=5, pady=5, sticky=tk.W)
            type_entry = ttk.Entry(self.keywd_frame, width=10)
            type_entry.grid(row=row, column=1, padx=5, pady=5)

            ttk.Label(self.keywd_frame, text="keywords").grid(row=row, column=2, padx=5, pady=5, sticky=tk.W)
            number_entry = ttk.Entry(self.keywd_frame, width=10)
            number_entry.grid(row=row, column=3, padx=5, pady=5)

            self.keywd_entries.append((type_entry, number_entry))

    def show_result(self, input_text):
        result_window = tk.Toplevel(self.root)
        result_window.title("Generated Input")
        result_window.geometry("400x300")

        text_box = tk.Text(result_window, wrap=tk.WORD)
        text_box.insert(tk.END, input_text)
        text_box.config(state=tk.DISABLED)
        text_box.pack(expand=True, fill=tk.BOTH, padx=10, pady=10)

        ttk.Button(result_window, text="Close", command=result_window.destroy).pack(pady=5)

#class Creat_xmi_input:
#    def __init__(self,root):
#
#    def generate_input(self):
#        # Collect standard inputs
#        self.fname=self.file_name.get()
#        for key, entry in self.entries.items():
#            try:
#                value = entry.get()
#                if value:
#                    self.ctrl_inputs[key]=value
#            except AttributeError:
#                raise TypeError(f"Expected tk.Entry, but got {type(entry)} for key {key}")
#        print('ctrl_inputs',self.ctrl_inputs)
#        return(self.ctrl_inputs)
#
#        # Collect fragment inputs if available
##        fragment_inputs = [self.description]
##        for pane in self.fragment_entries:
##            fragment_type = pane[0].get()
##            fragment_number = pane[1].get()
##            if fragment_type and fragment_number:
##                fragment_inputs.append((fragment_type, fragment_number))
##            else:
##                raise ValueError(f"Invalid pane in fragment_entries: {pane}")
##
#        # collect orbital inputs if available
#        if self.other_inp == True:
#            orbital_inputs = [self.description_orb]
#            for pane in self.orbital_entries:
#                orbital_number = pane.get()
#                if orbital_number:
#                    orbital_inputs.append(orbital_number)
#                else:
#                    raise ValueError(f"Invalid pane in orbital_entries: {pane}")
#
#        # collect spatial keyword inputs if available
#            keywd_inputs = [self.description_key]
#            for pane in self.keywd_entries:
#                keywd_type = pane[0].get()
#                keywd_number = pane[1].get()
#                if keywd_type and keywd_number:
#                    keywd_inputs.append((keywd_type, keywd_number))
#                else:
#                    raise ValueError(f"Invalid pane in keyword_entries: {pane}")
#
#            # Generate input string
#            self.input_text = f"{self.fname}\n"
#            self.input_text += "$ctrl\n"
#            self.input_text += " ".join(f"{key}={value}" for key, value in self.ctrl_inputs.items() if key != "fragments" and key != "orbital")
#            self.input_text += "\n$end"
#            self.input_text += "\n$frag\n"
#            self.input_text += f"{self.description}\n"
#            for frag in fragment_inputs[1:]:
#                self.input_text += f"{frag[0]} {frag[1]}\n"
#            self.input_text += "$end"
#            self.input_text += "\n$orb\n"
#            self.input_text += f"{self.description_orb}\n"
#            for orb in orbital_inputs[1:]:
#                self.input_text += f"{orb}\n"
#            self.input_text += "$end\n"
#            for keywd in keywd_inputs[1:]:
#                self.input_text += f"{keywd[0]} {keywd[1]}\n"
#
#            self.show_result(self.input_text)
#            
#            return(orbital_inputs, keywd_inputs)
#            
#    
#
#    def get_data(self):
#        input_text1=self.input_text
#        return input_text1
#    def get_file_name(self):
#        print('file_name',self.fname)
#        return self.fname

class analyse_inputs:
    def __init__(self, ctrl_data, fragment_data, orbital_data, spl_keywd_data):
        self.ctrl_data=ctrl_data
        self.fragment_data = fragment_data
        self.orbital_data = orbital_data 
        self.spl_keywd_inputs = spl_keywd_data

    def creat_input(self):
        self.keys_int = []       # List for ctrl keywd 
        self.values_int = []     # List for ctrl keywd's int input
        self.keys_str = []       # List for ctrl keywd
        self.values_str = []     # List for ctrl keywd's char input
        self.frag_type = []      # List for fragment type
        self.frag_atn = []       # List for associated atoms


        for i, (k, v) in enumerate(self.ctrl_data.items()):
            print('i', i, 'key', k, 'value', v)

            if i != 0 and i != 4:  # Check if the index is neither 0 nor 4
                self.keys_int.append(k)
                self.values_int.append(int(v))  # Convert to integer
            else:
                # Handle special cases for indices 0 and 4 (if any, as per your logic)
                self.keys_str.append(k)
                self.values_str.append(v)


        print('self.values_int',self.values_int)
        print('self.keys', self.keys_int)
        
    def run_Fortran(self):
        cisvb.cisvb_inp(self.keys_int, self.values_int, len(self.keys_int),
                        self.keys_str, self.values_str, len(self.keys_str) )
    
        if process.returncode != 0:
            print(f"Error running Fortran program: {stderr.decode()}")
        else:
            print(f"Fortran program executed successfully.\nOutput:\n{stdout.decode()}")

class Output:
    def __init__(self, root):
        # Initialize the root window
        self.root = root
        self.root.title("Output Value")

        # Create a frame for the widgets
        self.frame = ttk.Frame(root, padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Create a frame for the text widget and scrollbar
        text_frame = ttk.Frame(self.frame)
        text_frame.grid(row=1, column=0, pady=(0, 10), sticky=(tk.W, tk.E, tk.N, tk.S))

        # Add a vertical scrollbar
        self.scrollbar = ttk.Scrollbar(text_frame, orient=tk.VERTICAL)
        self.scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))

        # Add a Text widget to display the file content
        self.output_text = tk.Text(
            text_frame, wrap=tk.WORD, height=25, width=80, yscrollcommand=self.scrollbar.set
        )
        self.output_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

        # Configure the scrollbar to work with the Text widget
        self.scrollbar.config(command=self.output_text.yview)

    def display_file_content(self, output_file):
        """Displays the content of the specified file in the Text widget."""
        with open(output_file, 'r') as f:
            content = f.read()

        # Display the content in the Text widget
        self.output_text.delete(1.0, tk.END)  # Clear previous content
        self.output_text.insert(tk.END, content)

def open_second_root():
    # getting inputs and trasfer for analyse and creating input file
    ctrl_data,fragment_inputs, orbital_inputs, keywd_inputs=inputc.generate_input()
    print('ctrl_data',ctrl_data)
#    input_data=inputc.get_data()
#    file_name=inputc.get_file_name()
#    result=analyse_inputs(input_data)
    result=analyse_inputs(ctrl_data, fragment_inputs, orbital_inputs, keywd_inputs)
#    result.creat_input(file_name)
    result.creat_input()
    result.run_Fortran()     # instance to run fortran code

    # opening new Tkinter window and showing the results
    root1 = tk.Tk()
    output = Output(root1)
    output_file = "structure_set_1.dat"
    output.display_file_content(output_file)
    root1.mainloop()


def finish():
    sys.exit()

def create_orb(root):
    Orb_Input(root)

def create_keywd(root):
    Spl_Keywd(root)

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Input File Creator")
    root.geometry("500x400")
    inputc = Ctrl_Input(root)
    inputc.create_ctrl_pans()

    orbital_button = ttk.Button(root, text="Orbitals", command=create_orb(root))
    orbital_button.grid(row=1, column=0, padx=5, pady=10)

    keywd_button = ttk.Button(root, text="spatial keywords", command=create_keywd(root))
    keywd_button.grid(row=1, column=1, padx=5, pady=10)

    run_button = ttk.Button(root, text="RUN", command=open_second_root)
    run_button.grid(row=2, column=0, pady=10)

#    reset_button = ttk.Button(root, text="Reset", command=self.reset_fields)
#    reset_button.grid(row=12, column=2, padx=5, pady=10, sticky=tk.W)

    close_button = ttk.Button(root, text="FINISH", command=finish)
    close_button.grid(row=2, column=1, pady=10)

    root.mainloop()

#if __name__ == "__main__":
