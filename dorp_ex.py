import tkinter as tk
from tkinter import ttk

class DropdownMenuExample:
    def __init__(self, root):
        self.root = root
        self.root.title("Dropdown Entry Menu Example")
        self.root.geometry("400x300")
        self.root.configure(background="lightblue")

        # Label for dropdown
        label = ttk.Label(self.root, text="Select an Option:", style="Colour_Label.TLabel")
        label.grid(row=0, column=0, padx=10, pady=10, sticky=tk.W)

        # Dropdown menu options
        self.options = ["Option 1", "Option 2", "Option 3", "Option 4"]
        
        # StringVar to hold the selected value
        self.selected_option = tk.StringVar(value="Option 1")  # Default value

        # Create Combobox (dropdown)
        self.dropdown = ttk.Combobox(
            self.root,
            textvariable=self.selected_option,
            values=self.options,
            state="readonly"  # Makes the dropdown read-only
        )
        self.dropdown.grid(row=0, column=1, padx=10, pady=10)

        # Button to submit or confirm selection
        submit_button = ttk.Button(
            self.root,
            text="Submit",
            command=self.get_selected_option
        )
        submit_button.grid(row=1, column=0, columnspan=2, pady=20)

    def get_selected_option(self):
        # Function to handle selected option
        selected = self.selected_option.get()
        print(f"Selected Option: {selected}")
        tk.messagebox.showinfo("Selection", f"You selected: {selected}")

# Main application
if __name__ == "__main__":
    root = tk.Tk()
    app = DropdownMenuExample(root)
    root.mainloop()

