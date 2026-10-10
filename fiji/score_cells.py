#@ File (label="Scoring folder (<output>/scoring/<plate>_<marker>)", style="directory") folder
#@ String (label="Scorer initials", value="JY") scorer
"""
Blinded manual scoring of single-cell crops (Fiji, Jython; part of diffquant).
Make the crops first with:  diffquant score-sets config.yaml
Open this file in Fiji (File > Open, or drag onto Fiji) and press Run; choose
one scoring folder. Score each cell by eye as positive / negative for the
marker named in the window title - decide from the stain inside the outline,
not from neighbouring cells.

One cell at a time, fixed display range (never adjust B&C while scoring).
Keys (or buttons):
    P  positive      N  negative      U  unsure
    B  back (re-score previous)       O  outline on/off
Every answer is saved at once to calls_<scorer>.csv in the scoring folder;
re-running the script resumes at the first unscored cell.
"""
import csv
import os
import time

from ij import IJ, ImagePlus
from ij.gui import Overlay
from java.awt import BorderLayout, Font, GridLayout
from java.awt.event import KeyAdapter, KeyEvent, WindowAdapter
from javax.swing import JButton, JFrame, JLabel, JPanel, SwingUtilities

FOLDER = folder.getAbsolutePath()
CALLS = os.path.join(FOLDER, "calls_%s.csv" % scorer.strip())
KEYS = {KeyEvent.VK_P: "positive", KeyEvent.VK_N: "negative", KeyEvent.VK_U: "unsure"}

stack_imp = IJ.openImage(os.path.join(FOLDER, "crops.tif"))
stack = stack_imp.getStack()
n = stack.getSize()
all_rois = stack_imp.getOverlay()
vmin, vmax = 0.0, None
MARKER = os.path.basename(FOLDER).split("_")[-1]
with open(os.path.join(FOLDER, "display.txt")) as f:
    for line in f:
        if not line.strip() or line.startswith("#"):
            continue
        k, v = line.split()[:2]
        if k == "min":
            vmin = float(v)
        elif k == "max":
            vmax = float(v)
        elif k == "marker":
            MARKER = v

calls = {}
if os.path.exists(CALLS):
    with open(CALLS) as f:
        for row in csv.DictReader(f):
            calls[int(row["slice"])] = row


def save():
    with open(CALLS, "wb") as f:
        w = csv.writer(f)
        w.writerow(["slice", "call", "scorer", "seconds", "timestamp"])
        for s in sorted(calls):
            r = calls[s]
            w.writerow([s, r["call"], r["scorer"], r["seconds"], r["timestamp"]])


def first_unscored():
    for s in range(1, n + 1):
        if s not in calls:
            return s
    return None


state = {"slice": first_unscored() or n, "shown_at": time.time(), "outline": True}
view = ImagePlus("cell", stack.getProcessor(state["slice"]).duplicate())


def roi_for(s):
    for i in range(all_rois.size()):
        r = all_rois.get(i)
        if r.getPosition() == s:
            c = r.clone()
            c.setPosition(0)
            return c
    return None


def show(s):
    state["slice"] = s
    view.setProcessor(stack.getProcessor(s).duplicate())
    view.setDisplayRange(vmin, vmax)            # fixed for every cell
    ov = Overlay()
    r = roi_for(s)
    if r is not None and state["outline"]:
        ov.add(r)
    view.setOverlay(ov)
    view.updateAndDraw()
    done = len([k for k in calls if k <= n])
    status.setText("Cell %d / %d    scored %d    %s" % (
        s, n, done, ("previous answer: " + calls[s]["call"]) if s in calls else ""))
    state["shown_at"] = time.time()


def answer(call):
    s = state["slice"]
    calls[s] = {"call": call, "scorer": scorer,
                "seconds": "%.1f" % (time.time() - state["shown_at"]),
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S")}
    save()
    nxt = first_unscored()
    if nxt is None:
        status.setText("All %d cells scored - saved to %s" % (n, os.path.basename(CALLS)))
        IJ.showMessage("Scoring finished", "All %d cells scored.\nSaved: %s" % (n, CALLS))
        return
    show(nxt)


def back():
    if state["slice"] > 1:
        show(state["slice"] - 1)


def toggle_outline():
    state["outline"] = not state["outline"]
    show(state["slice"])


class Keys(KeyAdapter):
    def keyPressed(self, e):
        code = e.getKeyCode()
        if code in KEYS:
            answer(KEYS[code])
        elif code == KeyEvent.VK_B:
            back()
        elif code == KeyEvent.VK_O:
            toggle_outline()
        e.consume()


# --- windows ---------------------------------------------------------------
view.show()
win = view.getWindow()
canvas = win.getCanvas()
# Drop Fiji's single-key shortcuts on this window so P/N/U/B/O only score
for comp in (win, canvas):
    for kl in comp.getKeyListeners():
        comp.removeKeyListener(kl)
keys = Keys()
win.addKeyListener(keys)
canvas.addKeyListener(keys)
IJ.run(view, "Set... ", "zoom=200")
win = view.getWindow()

frame = JFrame("Score: %s  (%s)" % (os.path.basename(FOLDER), MARKER))
status = JLabel("", JLabel.CENTER)
status.setFont(Font("SansSerif", Font.PLAIN, 14))
buttons = JPanel(GridLayout(1, 5, 6, 6))
for label, fn in (("Positive (P)", lambda e: answer("positive")),
                  ("Negative (N)", lambda e: answer("negative")),
                  ("Unsure (U)", lambda e: answer("unsure")),
                  ("Back (B)", lambda e: back()),
                  ("Outline (O)", lambda e: toggle_outline())):
    b = JButton(label, actionPerformed=fn)
    b.setFocusable(False)                   # keep keys going to the frame
    buttons.add(b)
panel = JPanel(BorderLayout(6, 6))
panel.add(status, BorderLayout.NORTH)
panel.add(buttons, BorderLayout.CENTER)
frame.getContentPane().add(panel)
frame.addKeyListener(keys)
frame.pack()
frame.setLocation(win.getX(), win.getY() + win.getHeight() + 10)
frame.setVisible(True)

show(state["slice"])
if first_unscored() is None:
    status.setText("All cells already scored (%s). Use B to revise." % os.path.basename(CALLS))
win.toFront()
