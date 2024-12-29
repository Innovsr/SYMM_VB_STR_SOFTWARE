#################################################################################################################
# python GUI interface to create the input file to run the fortran symmetric structure generation symm_str code #
#################################################################################################################
    
import tkinter as tk
from tkinter import filedialog
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
        self.root = root
        self.input_text=''
        self.insert = False
        self.orbital_button = None
        self.structure_type_entry = False
        self.method_type_entry = False
        self.unit_type_entry = False
        self.readgeo = None
        self.file_path = None

        Pmw.initialise(root)
        self.balloon = Pmw.Balloon(root)
        style_colour_frame = ttk.Style()
        style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
        style_colour_label = ttk.Style()
        style_colour_label.configure("Colour_Label.TLabel", foreground="black", background="Ivory", relief="flat", font=("Arial", 14))
        style_colour_label1 = ttk.Style()
        style_colour_label1.configure("Colour_Label1.TLabel", foreground="black", background="lightblue", font=("Arial", 14))
        style = ttk.Style()
        style.configure("Custom.TRadiobutton", foreground="black", background="Lightblue", relief="flat", font=("Arial", 14))

        # Main Frame
        self.frame = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.frame1 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame1.grid(row=1, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.frame2 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame2.grid(row=2, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.frame3 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame3.grid(row=3, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.frame4 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
        self.frame4.grid(row=4, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        self.rect_border = ttk.Frame(self.frame2, relief="solid", borderwidth=2, padding=5)
        self.rect_border.grid(row=0, column=1)

        self.unit_type = tk.StringVar(value="None")
        # read prefered file name
        self.read_filename()
        self.create_geo_unit()

        # Input Fields
        self.entries = {}
        self.keywords = [" Number of Active Orbitals (nao) : ", " Number of Active Electrons (nae) : ", " Multiplicity of Active Part (nmul) : "]
        self.create_ctrl_pans()

        # read geometry from a file
        lebel = ttk.Label(self.frame1, text = "Brows to Upload Geometry File", style = "Colour_Label.TLabel")
        lebel.grid(row = 0, column=0, sticky=tk.W, padx=5, pady=5)
        button = ttk.Button(self.frame1, text = "Brows", command = self.read_geometry)
        button.grid(row = 0, column = 1, sticky = tk.W, padx = 5, pady = 5)

        # upload geometry manually
        lebel = ttk.Label(self.frame1, text = "Insert Geometry Manually", style = "Colour_Label.TLabel")
        lebel.grid(row = 2, column=0, sticky=tk.W, padx=5, pady=5)
        button = ttk.Button(self.frame1, text = "Geometry", command = self.insert_geo_manually)
        button.grid(row = 2, column = 1, sticky = tk.W, padx = 5, pady = 5)

        insert_button = ttk.Button(self.frame4, text="Insert", command=self.validate_and_generate)
        insert_button.grid(row=4, column=1, padx=5, pady=10)
        self.balloon.bind(insert_button, "Press Insert after providing all control keywords")

        self.tip_label = ttk.Label(self.frame3, text="Select structure type: ", style = "Colour_Label.TLabel")
        self.tip_label.grid(row = 0, column = 0, pady=10)
        self.balloon.bind(self.tip_label,"Select any one from the three types. select 'Covalent' to generate covalent sets only.\
Select 'Ionic' to generate ionic structures sets only. Select both to generate both type of structure sets.")
        self.tip_method_label = ttk.Label(self.frame3, text="Select method type: ", style = "Colour_Label.TLabel")
        self.tip_method_label.grid(row = 1, column = 0, pady=10)
        self.balloon.bind(self.tip_method_label,"Select any one from the two method provided. select 'Rumer' to genarate a Rumer structure\
set, no need to provide any spatial keywords. For all other type of set please select Chem Inst.")

        # Variable to store the selected input
        self.str_type = tk.StringVar(value="None")
        self.method_type = tk.StringVar(value="None")

        # Create radio buttons
        self.str_tip_buttons()

    def create_geo_unit(self):
        units = ["Bohr", "Angs"]
        label = ttk.Label(self.frame1, text = "Unit of the Geometry Data", style = "Colour_Label.TLabel")
        label.grid(row = 3, column = 0, sticky = tk.W, padx = 10, pady = 10)
        for i, unit in enumerate(units, start=1):
            button = ttk.Radiobutton(
                self.frame1,
                text=unit,
                value=unit,
                variable=self.unit_type,
                command=self.update_geo_unit,
                style="Custom.TRadiobutton"
            )
            button.grid(row=3, column=i, padx=10, pady=10)

    def update_geo_unit(self):
        geo_unit = self.unit_type.get()
        if geo_unit:
            self.unit_type_entry = True
            print('geo_unit',geo_unit)
            return (geo_unit)

    def read_filename(self):
        lebel1 = ttk.Label(self.frame, text="Name of the output file (optional)", style = "Colour_Label.TLabel").grid(row = 0, column=0, \
                sticky=tk.W, padx=5, pady=5)
        self.file_name = ttk.Entry(self.frame, width=20)
        self.file_name.grid(row = 0, column=1, padx=5, pady=5)
        if self.file_name == None:
            self.file_name='CISVB_OUT.DAT'
        return (self.file_name)

    def read_geometry(self):
        # Open file dialog to select a file
        self.file_path = filedialog.askopenfilename(title="Select Geometry File")

        if self.file_path:  # Check if a file is selected
            readgeo = Read_Geo(self.file_path)  # Initialize the Read_Geo class
            readgeo.read_geometry()       # Call the read_geometry method

            # Display the selected file path in the UI
            ttk.Label(
                self.frame1,
                text=f"Selected: {file_path}",
                style="Colour_Label1.TLabel"
            ).grid(row=1, column=0, columnspan=2)

    def insert_geo_manually(self):
        """Allows manual insertion of geometry data via Read_Geo."""
        readgeo = Read_Geo(self.file_path)
        readgeo.insert_geo(self.root)  # Call insert_geo method


    def create_ctrl_pans(self):
        lebel2 = ttk.Label(self.rect_border, text="Enter Control (ctrl) Keywords", font=("Arial", 14))
        lebel2.grid(row = 2, column=0, columnspan=2, pady=5)

        for idx, key in enumerate(self.keywords):
            ttk.Label(self.frame4, text=key, style = "Colour_Label.TLabel").grid(row = idx + 1, column=0, sticky=tk.W, padx=5, pady=5)
            entry = ttk.Entry(self.frame4, width=20)
            entry.grid(row = idx + 1, column=1, padx=5, pady=5)
            self.entries[key] = entry

    def validate_and_generate(self):
        self.ctrl_inputs = []
        if self.structure_type_entry == False:
            messagebox.showerror("Structure Type Missing", "Please Select a Structure Type among Covalent, Ionic, and Both.")
            return
        if self.method_type_entry == False:
            messagebox.showerror("Method Type Missing", "Please Select a Method between Chemical Insight and Rumer.")
            return
        if self.unit_type_entry == False:
            messagebox.showerror("Unit Type Missing", "Please Select a Unit between Bohr and Angs (angstrom).")
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
            messagebox.showerror("Validation Error", "Please insert valid values before insert other data.")
            return
        nao, nae, nmul = self.ctrl_inputs
        # Collect all inputs
#        print("Control Inputs:", self.ctrl_inputs)
#        for key, entry in self.entries.items():
#            self.ctrl_inputs[key] = entry.get().strip()
        return(nao)

    def str_tip_buttons(self):
        # Tips to display
        tips = ["Covalent", "Ionic", "Both"]
        tip_method = ["Chem inst", "Rumer"]

        # Create and pack each radio button
        i=2
        for tip in tips:
            i=i+1
            button = ttk.Radiobutton(
                self.frame3,
                text=tip,
                value=tip,
                variable=self.str_type,
                command=self.update_str_type,
                style="Custom.TRadiobutton"
            )
            button.grid(row = 0, column = i, padx=10, pady=10)

        i=2
        for tip in tip_method:
            i=i+1
            button = ttk.Radiobutton(
                self.frame3,
                text=tip,
                value=tip,
                variable=self.method_type,
                command=self.update_method_type,
                style="Custom.TRadiobutton"
            )
            button.grid(row = 1, column = i, padx=10, pady=10)

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

####################################################################################
######## Reading geometry starts here :
####################################################################################

class Read_Geo:
    def __init__(self, file_path):
        self.file_path = file_path
        self.geo_frame = None
        style_colour_label = ttk.Style()
        style_colour_label.configure("Colour_Label.TLabel", foreground="black", background="Ivory", relief="flat", font=("Arial", 14))
        style_colour_frame = ttk.Style()
        style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
        self.atoms = []  # List to store atom data


        
    def Deconvert(self, value):
        # Check if the string contains 'D' for scientific notation
        if 'D' in value:
            # Replace 'D' with 'e' for Python to interpret it
            value = value.replace('D', 'e')
        
        return float(value)


    def read_geometry(self):
#        Browse a file and read its lines to extract atom data.
#        Handles two formats of geometry files:
#        - 4 columns: atom name, x, y, z coordinates
#        - 5 columns: atom name, atomic number, x, y, z coordinates

        try:
            with open(self.file_path, 'r') as file:
                lines = file.readlines()


            for line in lines:
                line = line.strip()
                if not line:  # Skip empty lines
                    continue

                columns = line.split()
                if len(columns) == 4:
                    # 4-column format: atom name, x, y, z
                    atom_name = columns[0].capitalize()  # Capitalize atom name
                    x = self.Deconvert(columns[1])
                    y = self.Deconvert(columns[2])
                    z = self.Deconvert(columns[3])
#                    x, y, z = map(float, columns[1:])
                    self.atoms.append({
                        "atom": atom_name,
                        "x": x,
                        "y": y,
                        "z": z
                    })
                elif len(columns) == 5:
                    # 5-column format: atom name, atomic number, x, y, z
                    atom_name = columns[0].capitalize()  # Capitalize atom name
                    atomic_number = self.Deconvert(columns[1])
#                    x, y, z = map(float, columns[2:])
                    x = self.Deconvert(columns[2])
                    y = self.Deconvert(columns[3])
                    z = self.Deconvert(columns[4])
                    self.atoms.append({
                        "atom": atom_name,
                        "atomic_number": atomic_number,
                        "x": x,
                        "y": y,
                        "z": z
                    })
                else:
                    messagebox.showerror("Unknown Format", "Geometry file format is unknown Please Check help")
                    continue

            return atoms

        except Exception as e:
            print(f"An error occurred: {e}")
            return None

    def insert_geo(self, root):
        geo_window = tk.Toplevel(root)
        geo_window.title("Geometry")
        geo_window.geometry("650x650")
        geo_window.configure(background="lightblue")
        entry_widgets = []


        frame1 = ttk.Frame(geo_window, style = "Colour_Frame.TFrame", padding = 10 )
        frame1.grid(row = 0, column = 0)
        frame3 = ttk.Frame(geo_window, style = "Colour_Frame.TFrame", padding = 10 )
        frame3.grid(row = 2, column = 0)
         # Scrollable container
        container = ttk.Frame(geo_window)
        container.grid(row=1, column=0, sticky=tk.W+tk.E)

        canvas = tk.Canvas(container, background="lightblue", height=450, width=600)
        scrollbar = ttk.Scrollbar(container, orient=tk.VERTICAL, command=canvas.yview)
        frame2 = ttk.Frame(canvas, style="Colour_Frame.TFrame", padding=10)

        # Configure scrollable area
        frame2.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=frame2, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        label = ttk.Label(frame1, text = "Please insert the number of atoms: ", style ="Colour_Label.TLabel" )
        label.grid(row = 0, column = 0, padx = 10, pady = 10, sticky = tk.W)

        atom_num_entry = ttk.Entry(frame1, width = 10)
        atom_num_entry.grid(row = 0, column = 1, padx = 10, pady = 10 )

        def create_pan():
            atom_num = int(atom_num_entry.get())


            for widget in frame2.winfo_children():
                widget.destroy()
            entry_widgets.clear()


            label1 = ttk.Label(frame2, text = " Atom ", background ="lightblue")
            label1.grid(row = 1, column = 0, padx = 10)
            label2 = ttk.Label(frame2, text = " X Coordinates ", background ="lightblue")
            label2.grid(row = 1, column = 1, padx = 10)
            label3 = ttk.Label(frame2, text = " Y Coordinates ", background ="lightblue")
            label3.grid(row = 1, column = 2, padx = 10)
            label4 = ttk.Label(frame2, text = " Z Coordinates ", background ="lightblue")
            label4.grid(row = 1, column = 3, padx = 10)


            for i in range (atom_num):
                Atom = ttk.Entry(frame2, width = 10)
                Atom.grid(row = 2+i, column = 0, padx = 10, pady = 10)

                X_coord = ttk.Entry(frame2, width = 15)
                X_coord.grid(row = 2+i, column = 1, padx = 10, pady = 10)

                Y_coord = ttk.Entry(frame2, width = 15)
                Y_coord.grid(row = 2+i, column = 2, padx = 10, pady = 10)

                Z_coord = ttk.Entry(frame2, width = 15)
                Z_coord.grid(row = 2+i, column = 3, padx = 10, pady = 10)

                entry_widgets.append((Atom, X_coord, Y_coord, Z_coord))

        def fetch_data():
            for atom, x_entry, y_entry, z_entry in entry_widgets:
                atom_name = atom.get().strip()
                if atom:  # Only process if the atom field is not empty
                    x = float(self.Deconvert(x_entry.get().strip()))
                    y = float(self.Deconvert(y_entry.get().strip()))
                    z = float(self.Deconvert(z_entry.get().strip()))

                    self.atoms.append({
                        "atom": atom_name,
                        "x": x,
                        "y": y,
                        "z": z
                    })

            print("Atom Data:", self.atoms)
            return self.atoms


        button = ttk.Button(frame1, text = "Enter", command = create_pan)
        button.grid(row = 0, column = 3, padx = 10, pady = 10, sticky = tk.W )

        button = ttk.Button(frame3, text = "Insert", command = fetch_data)
        button.grid(row = 0, column = 3, padx = 10, pady = 10, sticky = tk.W )
        button = ttk.Button(frame3, text = "Close", command = geo_window.destroy)
        button.grid(row = 1, column = 3, padx = 10, pady = 10, sticky = tk.W )




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
                
#        style_colour_label1 = ttk.Style()
#        style_colour_label1.configure("Colour_Label1.TLabel", foreground="black", background="lightblue", relief="flat", font=("Arial", 16))


                # Insert button to validate and save inputs
        insert_button = ttk.Button(
            self.orbital_frame, text="Insert", command=self.validate_and_store_orbital_data
        )
        insert_button.grid(row=3, column=0, columnspan=2, pady=10)
        close_button = ttk.Button(self.orbital_frame, text = "Close", command = self.orbital_frame.destroy)
        close_button.grid(row = 4, column = 0, pady = 10, columnspan=2 )

    def create_orbital_section(self):
        # Create a new frame for orbital inputs if it doesn't exist
        if self.orbital_frame is None:
            self.orbital_frame = tk.Toplevel(root, padx=10, pady=10)
            self.orbital_frame.title("orbital inputs")
            self.orbital_frame.geometry("560x560")
            self.orbital_frame.configure(background="lightblue")

        for widget in self.orbital_frame.winfo_children():
            widget.destroy()

        frame = ttk.Frame(self.orbital_frame, style = "Colour_Frame.TFrame")
        frame.grid(row = 0, column = 0)
        frame0 = ttk.Frame(self.orbital_frame, style = "Colour_Frame.TFrame")
        frame0.grid(row = 1, column = 0)

         # Scrollable container
        container = ttk.Frame(self.orbital_frame)
        container.grid(row=2, column=0, sticky=tk.W+tk.E)

        canvas = tk.Canvas(container, background="lightblue", height=370, width=510)
        scrollbar = ttk.Scrollbar(container, orient=tk.VERTICAL, command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)

        frame1 = ttk.Frame(canvas, style = "Colour_Frame.TFrame")
        canvas.create_window((0, 0), window=frame1, anchor="nw")

        # Configure scrollable area
        frame1.bind(
            "<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=frame1, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)

        canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)
        scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

        # Create input panes for each active orbitals
        ttk.Label(frame, text=f"Number of active orbitals: {self.num_orbital}", style = "Colour_Label.TLabel" ).grid(row=0, column=0, columnspan=3, pady=15)

        #creating the pane for getting orbital numbers
        ttk.Label(frame0, text=f"     Atom Number", style = "Colour_Label1.TLabel").grid(row=0, column=3, padx=3, sticky=tk.E)
        ttk.Label(frame0, text=f"Sig-Pi / Fragment", style = "Colour_Label1.TLabel").grid(row=0, column=4, padx=5,  sticky=tk.E)

#        ttk.Label(frame1, text=f"--------------------", style = "Colour_Label1.TLabel").grid(row=2, column=1, padx=5, pady=5, sticky=tk.W)
#        ttk.Label(frame1, text=f"--------------------", style = "Colour_Label1.TLabel").grid(row=2, column=2, padx=5, pady=5, sticky=tk.W)

        for i in range(self.num_orbital):
            ttk.Label(frame1, text=f"active orbital {i+1} ", style = "Colour_Label.TLabel" ).grid(row=i+3, column=0, padx=30, pady=10, sticky=tk.W)
#       creating the pane for getting associated atom number
            self.assoatm_entry = ttk.Entry(frame1, width=10)
            self.assoatm_entry.grid(row=i+3, column=1, padx=30, pady=10)
            self.atm_entry.append(self.assoatm_entry)

#       creating the pane for getting associated atom number
            self.assotyp_entry = ttk.Entry(frame1, width=10)
            self.assotyp_entry.grid(row=i+3, column=2, padx=30, pady=10, sticky = tk.E)
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

class Keywd_Input:
    def __init__(self, root):
        self.root = root
        self.keywd_window = None
        self.priority_window = None
        self.priority_str_window = None
        self.frame = None
        self.set_type = tk.StringVar(value="None")
        self.cheminst_type = tk.StringVar(value="None")
        self.set_order_type = tk.StringVar(value="None")
        self.set_type_entry = False
        self.cheminst_type_entry = False
        self.set_order_type_entry = False
        Pmw.initialise(self.root)
        self.balloon = Pmw.Balloon(self.root)
        self.priority_values = {
            "IAB": "None",
            "SBB": "None",
            "NAB": "None",
            "PBU": "None",
            "PRU": "None",
        }


    def create_keywd_pane(self):
        if self.keywd_window is None:
            self.keywd_window = tk.Toplevel(self.root, padx = 10, pady = 10)
            self.keywd_window.title("Spatial Keyword Inputs")
            self.keywd_window.geometry("730x600")
            self.keywd_window.configure( background = "lightblue")

            self.frame = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame.grid (row = 0, column = 0, sticky = tk.W)
            self.frame1 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame1.grid (row = 1, column = 0, sticky = tk.W)
            self.frame9 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame9.grid (row = 9, column = 0)
            set_label1 = ttk.Label(self.frame, text = "Type of Output Set:", style = "Colour_Label.TLabel")
            set_label1.grid(row = 0, column = 0, sticky = tk.W, padx = 10, pady = 10)
            set_label2 = ttk.Label(self.frame1, text = "Type of Cheminst Set:", style = "Colour_Label.TLabel")
            set_label2.grid(row = 0, column = 0, padx = 10, pady = 10)
            close_button = ttk.Button(self.frame9, text = "Close", command = self.keywd_window.destroy)
            close_button.grid(row = 0, column = 0, pady = 10, columnspan=2 )

            Cheminst_type = ["Symmetry", "Asymmetry", "Checksymm"]
            for i, chem in enumerate(Cheminst_type, start=1):
                button1 = ttk.Radiobutton(
                    self.frame,
                    text=chem,
                    value=chem,
                    variable=self.cheminst_type,
                    command=self.cheminst_type_read,
                    style="Custom.TRadiobutton",
                )
                button1.grid(row=0, column=i, padx=10, pady=10)

            sets_type = ["Single Set", "All Best Sets", "All Sets", "Eq Bond"]
            for i, sets in enumerate(sets_type, start=1):
                button2 = ttk.Radiobutton(
                    self.frame1,
                    text=sets,
                    value=sets,
                    variable=self.set_type,
                    command=self.set_type_read,
                    style="Custom.TRadiobutton",
                )
                button2.grid(row=0, column=i, padx=10, pady=10)

            self.frame2 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame2.grid (row = 3, column = 0, sticky = tk.W)
            prio_button = ttk.Button(self.frame2, text = "Select Chemical Qualities According to Priorities", command = self.create_priority_keywd )
            prio_button.grid(row = 0, column = 0, padx = 10, pady = 5)
            self.balloon.bind(prio_button, "The sequence of the default priority is IAB > NAB > SBB and \n"
                              "PBU & PRU are not taken in the calculation; To change this default priorities\n"
                              "please click the button and selevct the numbers for each ")

            self.frame3 = ttk.Frame(self.keywd_window, style ="Colour_Frame.TFrame" )
            self.frame3.grid (row = 4, column = 0, sticky = tk.W)
            prio_struc_button = ttk.Button(self.frame3, text = "Insert Priorities Structures", command = self.Insert_Priority_Str )
            prio_struc_button.grid(row = 0, column = 0, padx = 10, pady = 5)
            self.balloon.bind(prio_struc_button, "If you want some structures must be present in the set \n"
                                                 "you can insert them by clicking this Button;\n" 
                                                 "The provided structures must be linearly independent \n"
                                                 "otherwise they will not be present in the set together")
    

    def set_type_read(self):
        set_type = self.set_type.get()
        if set_type:
            self.set_type_entry = True
            print('set_type',set_type)
            if set_type == 'All Sets':
                label = ttk.Label(self.frame1, text = "Enter Maximum Number of Output File:", style = "Colour_Label.TLabel")
                label.grid(row=1,column = 0, padx=10, pady=10, columnspan = 2)
                entry = ttk.Entry(self.frame1, width = 10)
                entry.grid(row = 1, column= 2, padx = 10, pady = 10)
            return (set_type)

    def cheminst_type_read(self):
        cheminst_type = self.cheminst_type.get()
        if cheminst_type:
            self.cheminst_type_entry = True
            print('cheminst_type',cheminst_type)
        #    return (cheminst_type)
            if cheminst_type == 'Symmetry':
                set_order = ["Big-to-Small","Small-to-Big","Qulity-Arrange"]
                for i, seto in enumerate(set_order, start=1):
                    button3 = ttk.Radiobutton(
                        self.frame,
                        text=seto,
                        value=seto,
                        variable=self.set_order_type,
                        command=self.set_order_type_read,
                        style="Custom.TRadiobutton",
                    )
                    button3.grid(row=1, column=i, padx=10, pady=10)

    def set_order_type_read(self):
        set_order_type = self.set_order_type.get()
        if set_order_type:
            self.set_order_type_entry = True
            print('set_order_type',set_order_type)
            return (set_order_type)

    def create_priority_keywd(self):
        if self.priority_window is None:
            self.priority_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.priority_window.title("Chem_Inst Qualities")
            self.priority_window.geometry("300x300")
            self.priority_window.configure(background="lightblue")

            frame = ttk.Frame(self.priority_window, style="Colour_Frame.TFrame")
            frame.grid(row=0, column=0, sticky=tk.W)
            frame1 = ttk.Frame(self.priority_window, style="Colour_Frame.TFrame")
            frame1.grid(row=1, column=0, sticky=tk.W)
            frame2 = ttk.Frame(self.priority_window, style="Colour_Frame.TFrame")
            frame2.grid(row=2, column=0)

            label = ttk.Label(
                frame,
                text="  Select Numbers from the  \nDropdown Menu or put None",
                style="Colour_Label.TLabel",
            )
            label.grid(row=0, column=0, padx=10, pady=5)

            # Priority options for dropdown
            Priority_Options = ["1", "2", "3", "4", "5", "None"]

            # initialise dictionary
            self.comboboxes = {}


            # Create labels and dropdowns for each priority
            labels = {
                "IAB": "Intra Atomic Bond: Minimize these bonds.",
                "SBB": "Symmetry Breaking Bond: Minimize these bonds.",
                "NAB": "Neighboring Atom Bond: Maximize these bonds.",
                "PBU": "Predefined Bond by User: Maximize these bonds.",
                "PRU": "Predefined Radicals by User: Maximize these structures.",
            }

            for i, (label_text, tooltip) in enumerate(labels.items(), start=1):
                # Create label
                label = ttk.Label(frame1, text=f"{label_text}:", style="Colour_Label1.TLabel")
                label.grid(row=i, column=0, padx=10, pady=5, sticky=tk.W)
                self.balloon.bind(label, tooltip)

                # Create combobox
                self.combobox = ttk.Combobox(frame1, values=Priority_Options, state="readonly", width=10)
                self.combobox.grid(row=i, column=1, padx=10, pady=5, sticky=tk.W)

                # Default selection
                self.combobox.set(self.priority_values[label_text])

                # Store the combobox reference
                self.comboboxes[label_text] = self.combobox


            close_button = ttk.Button(frame2, text = "DONE", command = self.get_priority_values)
            close_button.grid(row = 0, column = 0, pady = 10, columnspan=2 )
#            self.priority_window = False

    def get_priority_values(self):
        """
        Retrieve selected values from all priority dropdowns.
        """
        for key, self.combobox in self.comboboxes.items():
            self.priority_values[key] = self.combobox.get()
        print("Selected Priority Values:", self.priority_values)
        self.priority_window.destroy()
        self.priority_window = None

    def Insert_Priority_Str(self):
        if self.priority_str_window is None:
            self.priority_str_window = tk.Toplevel(self.root, padx=10, pady=10)
            self.priority_str_window.title("Priority Structure")
            self.priority_str_window.geometry("600x500")
            self.priority_str_window.configure(background="lightblue")

            frame = ttk.Frame(self.priority_str_window, style="Colour_Frame.TFrame")
            frame.grid(row=0, column=0, sticky=tk.W)
#            frame1 = ttk.Frame(self.priority_str_window, style="Colour_Frame.TFrame")
#            frame1.grid(row=1, column=0, sticky=tk.W)

             # Scrollable frame for dynamic fields
            canvas = tk.Canvas(self.priority_str_window, background="lightblue", height=300, width=560)
            scrollbar = ttk.Scrollbar(self.priority_str_window, orient=tk.VERTICAL, command=canvas.yview)
            scrollable_frame = ttk.Frame(canvas, style = "Colour_Frame.TFrame")

            # Configure scroll region
            scrollable_frame.bind(
                "<Configure>",
                lambda e: canvas.configure(scrollregion=canvas.bbox("all")),
            )
            canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
            canvas.configure(yscrollcommand=scrollbar.set)

            # Place the scrollable frame and scrollbar
            canvas.grid(row=1, column=0, sticky=tk.W + tk.E)
            scrollbar.grid(row=1, column=1, sticky=tk.N + tk.S)

            frame1 = ttk.Frame(self.priority_str_window, style="Colour_Frame.TFrame")
            frame1.grid(row=2, column=0)

            label = ttk.Label(frame, text = "Number of Structures:", style = "Colour_Label.TLabel")
            label.grid(row = 0, column = 0 , padx = 10, pady = 10)
            str_entry = ttk.Entry(frame, width = 10 )
            str_entry.grid(row = 0, column = 1)

            Insert_button = ttk.Button(frame, text = "Insert", command =lambda:self.generate_priority_str_fields(scrollable_frame, str_entry))
            Insert_button.grid(row = 0, column = 2, padx = 10, pady = 10 )

            close_button = ttk.Button(frame1, text = "DONE", command = self.get_priority_str_data)
            close_button.grid(row = 0, column = 0, pady = 10, columnspan=2 )

    def generate_priority_str_fields(self, frame1, str_entry):
        """
        Generate labels and entry fields dynamically based on user input.
        """
        try:
            str_number = int(str_entry.get())  # Get the number of structures
        except ValueError:
            # Handle invalid input gracefully
            tk.messagebox.showerror("Invalid Input", "Please enter a valid number.")
            return

        # Clear previous fields
        for widget in frame1.winfo_children():
            widget.destroy()

        # Create dynamic labels and entry fields
        self.structure_entries = []  # Store entry widget references
        for i in range(str_number):
            str_label = ttk.Label(frame1, text=f"Structure {i+1}:", style="Colour_Label1.TLabel")
            str_label.grid(row=i, column=0, padx=10, pady=10)

            struc_entry = ttk.Entry(frame1, width=50)
            struc_entry.grid(row=i, column=1, padx=10, pady=10)

            self.structure_entries.append(struc_entry)  # Save reference

    def get_priority_str_data(self):
        """
        Retrieve data from the dynamically created structure entry fields.
        """
        data = [entry.get() for entry in self.structure_entries]
        print("Entered Structures:", data)
        self.priority_str_window.destroy()
        self.priority_str_window = None
        return data







def finish():
    sys.exit()
        
def create_orb(inputc, root):
    num_orbital = inputc.generate_ctrl_input()  # calling generate_input to get active orbital number
    inputo=Orb_Input(num_orbital, root)

def create_keywd(root):
    input_keywd = Keywd_Input(root)
    input_keywd.create_keywd_pane()

#def create_keywd(root):
#    Spl_Keywd(root)

if __name__ == "__main__":
    root = tk.Tk()
    root.title("Input File Creator")
    root.geometry("550x750")
    root.configure(background="lightblue")
    inputc = Ctrl_Input(root)
    inputc.create_ctrl_pans()

    style_colour_frame = ttk.Style()
    style_colour_frame.configure("Colour_Frame.TFrame", background="lightblue")
    frame1 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
    frame1.grid(row=5, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    frame2 = ttk.Frame(root, style = "Colour_Frame.TFrame", padding="10")
    frame2.grid(row=6, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

    orbital_button = ttk.Button(frame1, text="ORBITALS", command=lambda:create_orb(inputc, root))
    orbital_button.grid(row=0, column=0, padx=5, pady=10)
    orbital_button.config(state=tk.DISABLED)  # Initially disable the button
    inputc.orbital_button = orbital_button

    keywd_button = ttk.Button(frame1, text="KEYWDS", command=lambda:create_keywd(root))
    keywd_button.grid(row=0, column=1, padx=5, pady=10)
#    keywd_button.config(state=tk.DISABLED)  # Initially disable the button
#    inputc.keywd_button = keywd_button

    Run_button = ttk.Button(frame1, text="RUN", command=lambda:create_orb(inputc, root))
    Run_button.grid(row=1, column=0, padx=5, pady=10)
#    Run_button.config(state=tk.DISABLED)  # Initially disable the button
#    inputc.Run_button = keywd_button

    close_button = ttk.Button(frame2, text="FINISH", command=finish)
    close_button.grid(row=0, column=1, pady=10)

    root.mainloop()
