import javax.swing.JLabel;
import javax.swing.SwingUtilities;

public class HP15C {

    private void createAndShowGUI() {
        CalcFrame frame = new CalcFrame();
        frame.setVisible(true);
    }

    public static void main(String[] args) {
        final HP15C app = new HP15C();
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {
                app.createAndShowGUI();
            }
        });
    }
}
