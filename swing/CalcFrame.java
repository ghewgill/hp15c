import java.awt.Color;
import java.awt.Component;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Toolkit;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.StringSelection;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.ByteArrayOutputStream;
import java.io.FileReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import javax.swing.AbstractAction;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.KeyStroke;
import javax.swing.SwingConstants;
import javax.swing.Timer;
import org.mozilla.javascript.Context;
import org.mozilla.javascript.Function;
import org.mozilla.javascript.FunctionObject;
import org.mozilla.javascript.ScriptableObject;
import org.mozilla.javascript.Scriptable;

class CalcFrame extends JFrame {
    private Map<Character, ImageIcon> pixmaps = new HashMap<Character, ImageIcon>();
    private JLabel digit[];
    private JLabel decimal[];
    private JLabel neg;
    private JLabel user;
    private JLabel f;
    private JLabel g;
    private JLabel trigmode;
    private JLabel complex;
    private JLabel prgm;
    private JLabel[] helplabels = new JLabel[40*3];

    private Context cx;
    private Scriptable scope;

    class AlertFunction extends ScriptableObject implements Function {
        public Object call(Context cx, Scriptable scope, Scriptable thisObj, Object[] args) {
            String s = (String) FunctionObject.convertArg(cx, scope, args[0], FunctionObject.JAVA_STRING_TYPE);
            JOptionPane.showMessageDialog(CalcFrame.this, s);
            return null;
        }
        public Scriptable construct(Context cx, Scriptable scope, Object[] args) {
            return null;
        }
        public String getClassName() {
            return "alert";
        }
    }

    static class Timeout extends Timer {
        private static Context cx;
        private static Scriptable scope;

        private Scriptable func;

        Timeout(int delay, final Scriptable func, boolean repeats) {
            super(delay, new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    Function f = (Function) func;
                    f.call(cx, scope, null, new Object[] {});
                }
            });
            setRepeats(repeats);
        }

        static void init(Context cx, Scriptable scope) {
            Timeout.cx = cx;
            Timeout.scope = scope;
        }

        static Scriptable setInterval(Scriptable f, int ms) {
            Timeout t = new Timeout(ms, f, true);
            t.start();
            return (Scriptable) Context.javaToJS(t, scope);
        }

        static Scriptable setTimeout(Scriptable f, int ms) {
            Timeout t = new Timeout(ms, f, false);
            t.start();
            return (Scriptable) Context.javaToJS(t, scope);
        }

        static void clearInterval(Scriptable s) {
            Timeout t = (Timeout) Context.jsToJava(s, Timeout.class);
            t.stop();
        }

        static void clearTimeout(Scriptable s) {
            Timeout t = (Timeout) Context.jsToJava(s, Timeout.class);
            t.stop();
        }
    }

    public class Window {
        public void console(String s) {
            System.out.println(s);
        }
    }

    public class Console {
        public void log(String s) {
            //System.out.println(s);
        }
    }

    public class Display {
        public void clear_digit(int i) {
            digit[i].setVisible(false);
        }

        public void clear_digits() {
            for (int i = 0; i < 10; i++) {
                digit[i].setVisible(false);
                decimal[i].setVisible(false);
            }
            neg.setVisible(false);
        }

        public void clear_shift() {
            f.setVisible(false);
            g.setVisible(false);
        }

        public void set_comma(int i) {
            decimal[i].setIcon(pixmaps.get(','));
            decimal[i].setVisible(true);
        }

        public void set_complex(boolean on) {
            complex.setVisible(on);
        }

        public void set_decimal(int i) {
            decimal[i].setIcon(pixmaps.get('.'));
            decimal[i].setVisible(true);
        }

        public void set_digit(int i, char d) {
            digit[i].setIcon(pixmaps.get(d));
            digit[i].setVisible(true);
        }

        public void set_neg() {
            neg.setVisible(true);
        }

        public void set_prgm(boolean on) {
            prgm.setVisible(on);
        }

        public void set_shift(String mode) {
            if (mode.equals("f")) {
                f.setVisible(true);
            } else if (mode.equals("g")) {
                g.setVisible(true);
            }
        }

        public void set_trigmode(String mode) {
            if (mode == null) {
                trigmode.setVisible(false);
            } else {
                trigmode.setText(mode);
                trigmode.setVisible(true);
            }
        }

        public void set_user(boolean on) {
            user.setVisible(on);
        }
    }

    byte[] loadFile(InputStream in) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try {
            byte[] buf = new byte[16384];
            while (true) {
                int n = in.read(buf, 0, buf.length);
                if (n <= 0) {
                    break;
                }
                baos.write(buf, 0, n);
            }
        } catch (IOException e) {
            System.err.println(e);
        }
        return baos.toByteArray();
    }

    void jsCall(String name, Object... args) {
        Function fn = (Function) cx.evaluateString(scope, name, null, 1, null);
        fn.call(cx, scope, null, args);
    }

    public CalcFrame() {
        super("HP 15C");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        final JPanel outerpane = new JPanel();
        setContentPane(outerpane);
        outerpane.setLayout(null);

        final JPanel pane = new JPanel();
        outerpane.add(pane);
        pane.setLayout(null);

        int keymask = ActionEvent.CTRL_MASK;
        if (System.getProperty("os.name").equals("Mac OS X")) {
            keymask = ActionEvent.META_MASK;
        }

        JMenuBar menubar = new JMenuBar();
        JMenu editmenu = new JMenu("Edit");
        menubar.add(editmenu);
        JMenuItem copyitem = new JMenuItem("Copy");
        copyitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_C, keymask));
        copyitem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                String s = cx.evaluateString(scope, "Stack[0]", null, 1, null).toString();
                Toolkit.getDefaultToolkit().getSystemClipboard().setContents(new StringSelection(s), null);
            }
        });
        editmenu.add(copyitem);
        JMenuItem pasteitem = new JMenuItem("Paste");
        pasteitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_V, keymask));
        pasteitem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                try {
                    Object s = Toolkit.getDefaultToolkit().getSystemClipboard().getData(DataFlavor.stringFlavor);
                    jsCall("paste", s.toString());
                } catch (UnsupportedFlavorException x) {
                } catch (IOException x) {
                }
            }
        });
        editmenu.add(pasteitem);
        JMenu viewmenu = new JMenu("View");
        menubar.add(viewmenu);
        final JCheckBoxMenuItem keysitem = new JCheckBoxMenuItem("Full Keyboard");
        keysitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_K, keymask));
        keysitem.setSelected(true);
        viewmenu.add(keysitem);
        JMenu testmenu = new JMenu("Test");
        menubar.add(testmenu);
        JMenuItem testitem = new JMenuItem("Test");
        testitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T, keymask));
        testitem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                jsCall("start_tests");
            }
        });
        testmenu.add(testitem);
        setJMenuBar(menubar);

        final ImageIcon face = new ImageIcon(loadFile(getClass().getResourceAsStream("/15.jpg")));
        JLabel facelabel = new JLabel(face);
        pane.add(facelabel);
        facelabel.setBounds(0, 0, face.getIconWidth(), face.getIconHeight());
        pane.setBounds(0, 0, face.getIconWidth(), face.getIconHeight());

        keysitem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (keysitem.getState()) {
                    pane.setLocation(0, 0);
                    outerpane.setPreferredSize(new Dimension(face.getIconWidth(), face.getIconHeight()));
                    pack();
                } else {
                    pane.setLocation(-125, -50);
                    outerpane.setPreferredSize(new Dimension(350, 80));
                    pack();
                }
            }
        });

        final String pixmap_chars = "0123456789-ABCDEoru";
        for (int i = 0; i < pixmap_chars.length(); i++) {
            char c = pixmap_chars.charAt(i);
            pixmaps.put(c, new ImageIcon(loadFile(getClass().getResourceAsStream("/"+c+".png"))));
        }
        pixmaps.put('.', new ImageIcon(loadFile(getClass().getResourceAsStream("/"+"decimal.png"))));
        pixmaps.put(',', new ImageIcon(loadFile(getClass().getResourceAsStream("/"+"comma.png"))));

        digit = new JLabel[10];
        decimal = new JLabel[10];
        for (int i = 0; i < 10; i++) {
            digit[i] = new JLabel();
            pane.add(digit[i], 0);
            digit[i].setBounds(175 + i * 27, 67, 20, 30);
            decimal[i] = new JLabel();
            pane.add(decimal[i], 0);
            decimal[i].setBounds(194 + i * 27, 91, 6, 10);
            decimal[i].setVerticalAlignment(SwingConstants.TOP);
        }

        neg = new JLabel(new ImageIcon(loadFile(getClass().getResourceAsStream("/neg.png"))));
        pane.add(neg, 0);
        neg.setBounds(158, 80, 12, 3);

        Font font = new Font("sans-serif", Font.PLAIN, 10);

        user = new JLabel("USER");
        user.setFont(font);
        pane.add(user, 0);
        user.setBounds(190, 100, 30, 10);

        f = new JLabel("f");
        f.setFont(font);
        pane.add(f, 0);
        f.setBounds(230, 100, 30, 10);

        g = new JLabel("g");
        g.setFont(font);
        pane.add(g, 0);
        g.setBounds(250, 100, 30, 10);

        trigmode = new JLabel();
        trigmode.setFont(font);
        pane.add(trigmode, 0);
        trigmode.setBounds(300, 100, 30, 10);

        complex = new JLabel("C");
        complex.setFont(font);
        pane.add(complex, 0);
        complex.setBounds(390, 100, 30, 10);

        prgm = new JLabel("PRGM");
        prgm.setFont(font);
        pane.add(prgm, 0);
        prgm.setBounds(410, 100, 30, 10);

        cx = Context.enter();
        scope = cx.initStandardObjects();

        Timeout.init(cx, scope);
        try {
            ScriptableObject.putProperty(scope, "setTimeout", new FunctionObject("setTimeout", Timeout.class.getDeclaredMethod("setTimeout", new Class[] {Scriptable.class, int.class}), scope));
            ScriptableObject.putProperty(scope, "clearTimeout", new FunctionObject("clearTimeout", Timeout.class.getDeclaredMethod("clearTimeout", new Class[] {Scriptable.class}), scope));
            ScriptableObject.putProperty(scope, "setInterval", new FunctionObject("setInterval", Timeout.class.getDeclaredMethod("setInterval", new Class[] {Scriptable.class, int.class}), scope));
            ScriptableObject.putProperty(scope, "clearInterval", new FunctionObject("clearInterval", Timeout.class.getDeclaredMethod("clearInterval", new Class[] {Scriptable.class}), scope));
            ScriptableObject.putProperty(scope, "alert", new AlertFunction());
        } catch (NoSuchMethodException e) {
            System.err.println(e);
        }

        try {
            String[] files = {
                "/sprintf-0.6.js",
                "/matrix.js",
                "/hp15c.js",
                "/test.js",
            };
            for (String fn : files) {
                cx.evaluateReader(scope, new InputStreamReader(getClass().getResourceAsStream(fn)), fn, 1, null);
            }
        } catch (IOException e) {
        }

        for (int r = 0; r < 4; r++) {
            for (int c = 0; c < 10; c++) {
                int i = r * 10 + c;
                int th = 34;
                if (c == 5 && r >= 2) {
                    if (r == 2) {
                        th = 99;
                    } else {
                        continue;
                    }
                }
                final int bx = 81 + c * 57;
                final int by = 169 + r * 65;
                final int h = th;
                ImageIcon bicon = new ImageIcon(face.getImage()) {
                    public void paintIcon(Component c, Graphics g, int x, int y) {
                        g.drawImage(getImage(), 0, 0, 39, h, bx, by, bx+39, by+h, null);
                    }
                };
                JButton b = new JButton(bicon);
                b.setBorderPainted(false);
                b.setPressedIcon(new ImageIcon(face.getImage()) {
                    public void paintIcon(Component c, Graphics g, int x, int y) {
                        g.drawImage(getImage(), 0, 0, 39, h, bx, by+1, bx+39, by+1+h, null);
                    }
                });
                final String key = cx.evaluateString(scope, String.format("KeyTable[%1$d][%2$d]", r, c), null, 1, null).toString();
                b.addActionListener(new ActionListener() {
                    public void actionPerformed(ActionEvent e) {
                        jsCall("key", key);
                    }
                });
                b.setBounds(bx, by, 39, h);
                pane.add(b, 0);
                if (!(r == 3 && (c == 0 || c == 5))) {
                    String hk = key;
                    if (hk == "\b") {
                        hk = "\u2190";
                    } else if (hk == "\r") {
                        hk = "\u21b2";
                    }
                    JLabel help = new JLabel(hk);
                    help.setBounds(70 + 57 * c, 167 + 65 * r, 16, 16);
                    help.setOpaque(true);
                    help.setBackground(Color.yellow);
                    help.setHorizontalAlignment(SwingConstants.CENTER);
                    help.setFont(new Font("Courier", 0, 14));
                    help.setVisible(false);
                    pane.add(help, 0);
                    helplabels[i] = help;
                }
            }
        }
        Color goldenrod = new Color(218, 165, 32);
        Color lightblue = new Color(173, 216, 230);
        int i = 0;
        while (true) {
            Object inf = cx.evaluateString(scope, String.format("ExtraKeyTable[%1$d]", i), null, 1, null);
            if (inf == Context.getUndefinedValue()) {
                break;
            }
            Scriptable info = (Scriptable) inf;
            int r = ((Number) info.get(0, info)).intValue();
            int c = ((Number) info.get(1, info)).intValue();
            int f = ((Number) info.get(2, info)).intValue();
            String s = info.get(3, info).toString();
            int top = 167 + 65*r + 20*f;
            int left = 70 + 57*c;
            JLabel help = new JLabel(s);
            help.setBounds(left, top, 16, 16);
            help.setOpaque(true);
            help.setBackground(f == 1 ? lightblue : goldenrod);
            help.setHorizontalAlignment(SwingConstants.CENTER);
            help.setFont(new Font("Courier", 0, 14));
            help.setVisible(false);
            pane.add(help, 0);
            helplabels[40+r*10+c+(f>0 ? 40 : 0)] = help;
            i++;
        }

        Display display = new Display();

        display.clear_digits();
        display.set_user(false);
        display.clear_shift();
        display.set_trigmode(null);

        ScriptableObject.putProperty(scope, "Display", Context.javaToJS(new Display(), scope));
        ScriptableObject.putProperty(scope, "window", Context.javaToJS(new Window(), scope));
        ScriptableObject.putProperty(scope, "console", Context.javaToJS(new Console(), scope));

        jsCall("init");

        pane.getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(KeyStroke.getKeyStroke('h'), "help");
        pane.getActionMap().put("help", new AbstractAction() {
            public void actionPerformed(ActionEvent e) {
                for (JLabel hl : helplabels) {
                    if (hl != null) {
                        hl.setVisible(!hl.isVisible());
                    }
                }
            }
        });
        Scriptable chartable = (Scriptable) cx.evaluateString(scope, "CharTable", null, 1, null);
        for (Object id : chartable.getIds()) {
            final String k = id.toString();
            pane.getInputMap(JComponent.WHEN_ANCESTOR_OF_FOCUSED_COMPONENT).put(KeyStroke.getKeyStroke(k.charAt(0)), k);
            pane.getActionMap().put(k, new AbstractAction() {
                public void actionPerformed(ActionEvent e) {
                    jsCall("key", k);
                }
            });
        }

        //Display the window.
        outerpane.setPreferredSize(new Dimension(face.getIconWidth(), face.getIconHeight()));
        setResizable(false);
        pack();
    }
}
