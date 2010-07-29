import java.awt.Container;
import java.awt.Font;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyAdapter;
import java.awt.event.KeyEvent;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import javax.swing.ImageIcon;
import javax.swing.JCheckBoxMenuItem;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
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

        Timeout(int delay, final Scriptable func) {
            super(delay, new ActionListener() {
                public void actionPerformed(ActionEvent e) {
                    Function f = (Function) func;
                    f.call(cx, scope, null, new Object[] {});
                }
            });
            setRepeats(false);
        }

        static void init(Context cx, Scriptable scope) {
            Timeout.cx = cx;
            Timeout.scope = scope;
        }

        static Scriptable setTimeout(Scriptable f, int ms) {
            Timeout t = new Timeout(ms, f);
            t.start();
            return (Scriptable) Context.javaToJS(t, scope);
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

    byte[] loadFile(String fn) {
        ByteArrayOutputStream baos = new ByteArrayOutputStream();
        try {
            FileInputStream f = new FileInputStream(fn);
            byte[] buf = new byte[16384];
            while (true) {
                int n = f.read(buf, 0, buf.length);
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

    public CalcFrame() {
        super("HP 15C");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        Container pane = getContentPane();
        pane.setLayout(null);

        JMenuBar menubar = new JMenuBar();
        JMenu editmenu = new JMenu("Edit");
        menubar.add(editmenu);
        JMenuItem copyitem = new JMenuItem("Copy");
        copyitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_COPY, 0));
        editmenu.add(copyitem);
        JMenuItem pasteitem = new JMenuItem("Paste");
        pasteitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_PASTE, 0));
        editmenu.add(pasteitem);
        JMenu viewmenu = new JMenu("View");
        menubar.add(viewmenu);
        JCheckBoxMenuItem keysitem = new JCheckBoxMenuItem("Full Keyboard");
        keysitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_K, ActionEvent.CTRL_MASK));
        keysitem.setSelected(true);
        viewmenu.add(keysitem);
        JMenu testmenu = new JMenu("Test");
        menubar.add(testmenu);
        JMenuItem testitem = new JMenuItem("Test");
        testitem.setAccelerator(KeyStroke.getKeyStroke(KeyEvent.VK_T, ActionEvent.CTRL_MASK));
        testitem.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent e) {
                Object r = cx.evaluateString(scope, "start_tests()", null, 1, null);
            }
        });
        testmenu.add(testitem);
        setJMenuBar(menubar);

        ImageIcon face = new ImageIcon(loadFile("../15.jpg"));
        JLabel facelabel = new JLabel(face);
        pane.add(facelabel);
        facelabel.setBounds(0, 0, face.getIconWidth(), face.getIconHeight());

        final String pixmap_chars = "0123456789-ABCDEoru";
        for (int i = 0; i < pixmap_chars.length(); i++) {
            char c = pixmap_chars.charAt(i);
            pixmaps.put(c, new ImageIcon(loadFile("../"+c+".png")));
        }
        pixmaps.put('.', new ImageIcon(loadFile("../"+"decimal.png")));
        pixmaps.put(',', new ImageIcon(loadFile("../"+"comma.png")));

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

        neg = new JLabel(new ImageIcon(loadFile("../neg.png")));
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
            ScriptableObject.putProperty(scope, "alert", new AlertFunction());
        } catch (NoSuchMethodException e) {
            System.err.println(e);
        }

        try {
            String[] files = {
                "../sprintf-0.6.js",
                "../jsmat/matrix.js",
                "../hp15c.js",
                "../test.js",
            };
            for (String fn : files) {
                cx.evaluateReader(scope, new FileReader(fn), fn, 1, null);
            }
        } catch (IOException e) {
        }
        Display display = new Display();

        display.clear_digits();
        display.set_user(false);
        display.clear_shift();
        display.set_trigmode(null);

        Object d = Context.javaToJS(new Display(), scope);
        ScriptableObject.putProperty(scope, "Display", d);
        ScriptableObject.putProperty(scope, "window", Context.javaToJS(new Window(), scope));
        ScriptableObject.putProperty(scope, "console", Context.javaToJS(new Console(), scope));

        Object r = cx.evaluateString(scope, "init()", null, 1, null);

        addKeyListener(new KeyAdapter() {
            public void keyTyped(KeyEvent e) {
                String k = Character.toString(e.getKeyChar());
                switch (e.getKeyChar()) {
                    case '\n':
                        k = "\\r";
                        break;
                    case 0x14: // ^T
                        return;
                }
                Object r = cx.evaluateString(scope, "key('"+k+"')", null, 1, null);
            }
        });

        //Display the window.
        //frame.pack();
        setSize(face.getIconWidth(), face.getIconHeight());
    }
}
