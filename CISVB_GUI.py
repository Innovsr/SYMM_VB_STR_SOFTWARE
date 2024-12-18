#################################################################################################################
# python GUI interface to create the input file to run the fortran symmetric structure generation symm_str code #
#################################################################################################################
    
import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
#import tkinter.messagebox as msg
import re
import os
#import program_main
import subprocess
import sys
import cisvb
import numpy as np
import Pmw


class Ctrl_Input:
    def __init__(self, root):
# initialise disctionary to store control data
        self.ctrl_inputs={}
        self.input_text=''
        self.insert = False
        self.orbital_button = None
        self.structure_type_entry = False
        self.method_type_entry = False
        Pmw.initialise(root)
        self.balloon = Pmw.Balloon(root)
        style_colour_frame = ttk.Style()
        style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
        style_colour_label = ttk.Style()
        style_colour_label.configure("Colour_Label.TLabel", foreground="white", background="blue", font=("Arial", 14))

        # Main Frame
        self.frame = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.rect_border = ttk.Frame(self.frame, relief="solid", borderwidth=2, padding=5)
        self.rect_border.grid(row=2, column=1)

        # Input Fields
        self.entries = {}
        self.keywords = ["Number of Active Orbitals (nao) :          ", "Number of Active Orbitals (nae) :          ", "Multiplicity of Active Part (nmul) :          "]
        self.create_ctrl_pans()

        insert_button = ttk.Button(self.frame, text="Insert", command=self.validate_and_generate)
        insert_button.grid(row=9, column=1, padx=5, pady=10)
        self.balloon.bind(insert_button, "Press Insert after providing all control keywords")

        self.tip_label = ttk.Label(self.frame, text="Select structure type: ", style = "Colour_Label.TLabel")
        self.tip_label.grid(row = 3, column = 0, pady=10)
        self.balloon.bind(self.tip_label,"Select any one from the three types. select 'Covalent' to generate covalent sets only.\
Select 'Ionic' to generate ionic structures sets only. Select both to generate both type of structure sets.")
        self.tip_method_label = ttk.Label(self.frame, text="Select method type: ", style = "Colour_Label.TLabel")
        self.tip_method_label.grid(row = 4, column = 0, pady=10)
        self.balloon.bind(self.tip_method_label,"Select any one from the two method provided. select 'Rumer' to genarate a Rumer structure\
set, no need to provide any spatial keywords. For all other type of set please select Chem Inst.")

        # Variable to store the selected tip
        self.str_type = tk.StringVar(value="None")
        self.method_type = tk.StringVar(value="None")

        # Create radio buttons
        self.str_tip_buttons()


    def create_ctrl_pans(self):
        lebel1 = ttk.Label(self.frame, text="Name of the File (optional)", font=("Arial", 12)).grid(row =  1, column=0, sticky=tk.W, padx=5, pady=5)
        lebel2 = ttk.Label(self.rect_border, text="Enter Control (ctrl) Keywords", font=("Arial", 14))
        lebel2.grid(row = 2, column=0, columnspan=2, pady=5)

        self.file_name = ttk.Entry(self.frame, width=20)
        self.file_name.grid(row = 1, column=1, padx=5, pady=5)
        #        self.entries["file_name"]=entry
        for idx, key in enumerate(self.keywords):
            ttk.Label(self.frame, text=key, font=("Arial", 12)).grid(row = idx + 5, column=0, sticky=tk.W, padx=5, pady=5)
            entry = ttk.Entry(self.frame, width=20)
            entry.grid(row = idx + 5, column=1, padx=5, pady=5)
            self.entries[key] = entry

    def validate_and_generate(self):
        self.ctrl_inputs = []
        if self.structure_type_entry == False:
            messagebox.showerror("Structure Type Missing", "Please Select a Structure Type among Covalent, Ionic, and Both.")
            return
        if self.method_type_entry == False:
            messagebox.showerror("Method Type Missing", "Please Select a Method between Chemical Insight and Rumer.")
            return
        else:
            # Loop through each input field to validate
            for key, entry in self.entries.items():
                value = entry.get().strip()
                if not value:  # If any field is empty
                    messagebox.showerror("Validation Error", f"'{key}' is empty. Please fill in all fields.")
                    entry.configure(background="pink")  # Highlight the empty field
                    entry.focus_set()  # Focus the empty field
                    self.insert = False
                    return False
                else:
                    self.ctrl_inputs.append(value)

            for entry in self.entries.values():
                entry.configure(background="white")  # Reset background color
            self.insert = True
            if self.orbital_button:  # Enable the Orbital button
                nao, nae, nmul = self.ctrl_inputs
                print('nao, nae, nmul',nao, nae, nmul)
                self.orbital_button.config(state=tk.NORMAL)
            return True  # Validation successful

    def generate_ctrl_input(self):
        if not self.insert:
            messagebox.showerror("Validation Error", "Please insert valid values before inserting control input.")
            return
        # Collect all inputs
        for key, entry in self.entries.items():
            self.ctrl_inputs[key] = entry.get().strip()
        print("Control Inputs:", self.ctrl_inputs)
        return(self.ctrl_inputs)

    def str_tip_buttons(self):
        # Tips to display
        tips = ["Covalent", "Ionic", "Both"]
        tip_method = ["Chem inst", "Rumer"]
        style = ttk.Style()
        style.configure("Custom.TRadiobutton", font=("Arial", 12))

        # Create and pack each radio button
        i=0
        for tip in tips:
            i=i+1
            button = ttk.Radiobutton(
                self.frame,
                text=tip,
                value=tip,
                variable=self.str_type,
                command=self.update_str_type,
                style="Custom.TRadiobutton"
            )
            button.grid(row = 3, column = i, padx=10, pady=5)

        i=0
        for tip in tip_method:
            i=i+1
            button = ttk.Radiobutton(
                self.frame,
                text=tip,
                value=tip,
                variable=self.method_type,
                command=self.update_method_type,
                style="Custom.TRadiobutton"
            )
            button.grid(row = 4, column = i, padx=10, pady=5)

    def update_str_type(self):
        structure_type = self.str_type.get()
        if structure_type:
            self.structure_type_entry = True
            print('structure_type',structure_type)
            return (structure_type)

    def update_method_type(self):
        method_type = self.method_type.get()
        if method_type:
            self.method_type_entry = True
            print('method_type',method_type)
            return (method_type)





###################################################################################
########## orbital section starts here :
###################################################################################

class Orb_Input:
    def __init__(self,num_orbital, root):
        self.num_orbital = int(num_orbital)

        self.orbital_frame = None
        self.atm_entry = []
        self.typ_entry = []
#       self.assoatm_entries = []
#       self.assotyp_entries = []
        self.orbital_data = []
        self.description_orb = ""
        self.create_orbital_section()
                
                # Insert button to validate and save inputs
        insert_button = ttk.Button(
            self.orbital_frame, text="Insert", command=self.validate_and_store_orbital_data
        )
        insert_button.grid(row=self.num_orbital + 4, column=1, columnspan=2, pady=10)

    def create_orbital_section(self):
        # Create a new frame for orbital inputs if it doesn't exist
        if self.orbital_frame is None:
            self.orbital_frame = tk.Toplevel(root, padx=10, pady=10)
            self.orbital_frame.title("orbital inputs")
            self.orbital_frame.geometry("600x600")

        for widget in self.orbital_frame.winfo_children():
            widget.destroy()


        # Create input panes for each active orbitals
        ttk.Label(self.orbital_frame, text=f"Number of active orbitals: {self.num_orbital}").grid(
            row=0, column=0, columnspan=2, pady=5
        )

        #creating the pane for getting orbital numbers
        ttk.Label(self.orbital_frame, text=f"associated atom number").grid(row=2, column=3, padx=5, pady=5, sticky=tk.W)
        ttk.Label(self.orbital_frame, text=f"Sig-Pi/fragment").grid(row=2, column=4, padx=5, pady=5, sticky=tk.W)

        for i in range(self.num_orbital):
            ttk.Label(self.orbital_frame, text=f"active orbital {i+1} ").grid(row=i+3, column=1, padx=5, pady=5, sticky=tk.W)
#       creating the pane for getting associated atom number
            self.assoatm_entry = ttk.Entry(self.orbital_frame, width=10)
            self.assoatm_entry.grid(row=i+3, column=3, padx=5, pady=5)
            self.atm_entry.append(self.assoatm_entry)

#       creating the pane for getting associated atom number
            self.assotyp_entry = ttk.Entry(self.orbital_frame, width=10)
            self.assotyp_entry.grid(row=i+3, column=4, padx=5, pady=5)
            self.typ_entry.append(self.assotyp_entry)


    def validate_and_store_orbital_data(self):
        """Validate orbital inputs and store them in a list."""
        self.orbital_data.clear()  # Clear previous data

        for i in range(self.num_orbital):
            atom_number = self.atm_entry[i].get().strip()
            orbital_type = self.typ_entry[i].get().strip()

            # Validate inputs
            if not atom_number or not orbital_type:
                messagebox.showerror(
                    "Validation Error",
                    f"Please fill in all fields for orbital {i + 1}.",
                )
                self.atm_entry[i].configure(background="pink")
                self.typ_entry[i].configure(background="pink")
                return

            try:
                atom_number = int(atom_number)  # Ensure atom number is an integer
            except ValueError:
                messagebox.showerror(
                    "Input Error", f"Atom number for orbital {i + 1} must be an integer."
                )
                self.atm_entry[i].configure(background="pink")
                return

            # Reset background color after validation
            self.atm_entry[i].configure(background="white")
            self.typ_entry[i].configure(background="white")

            # Store the data
            self.orbital_data.append({
                "atom_number": atom_number,
                "orbital_type": orbital_type
            })

        messagebox.showinfo("Success", "Orbital inputs validated and stored successfully.")



        # collect orbital inputs if available
#    def generate_orb_input(self):
#        for pane in self.atm_entry:
#            input_orb = pane.get()
#            if input_orb:
#                self.assoatm_entries.append(input_orb)
#            else:
#                raise ValueError(f"Invalid pane in associated atom number entries: {pane}")
#        for pane in self.typ_entry:
#            input_orb = pane.get()
#            if input_orb:
#                self.assotyp_entries.append(input_orb)
#            else:
#                raise ValueError(f"Invalid pane in sig-pi/fragment entries: {pane}")
#        print ('self.assoatm_entries, self_assotyp_entries',self.assoatm_entries, self_assotyp_entries)
#        return (self.assoatm_entries, self_assotyp_entries)



def finish():
    sys.exit()
        
def create_orb(inputc, root):
    ctrl_inputs = []
    ctrl_inputs = inputc.generate_ctrl_input()  # calling generate_input to get ctrl values
    num_orbital = ctrl_inputs["nao"]  # read the number of active orbitals
    inputo=Orb_Input(num_orbital, root)

#def create_keywd(root):
#    Spl_Keywd(root)

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Input File Creator")
    root.geometry("800x550")
    root.configure(background="lightblue")
    inputc = Ctrl_Input(root)
    inputc.create_ctrl_pans()

    orbital_button = ttk.Button(root, text="ORBITALS", command=lambda:create_orb(inputc, root))
    orbital_button.grid(row=11, column=0, padx=5, pady=10)
    orbital_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.orbital_button = orbital_button

    close_button = ttk.Button(root, text="FINISH", command=finish)
    close_button.grid(row=13, column=0, pady=10)

    root.mainloop()
